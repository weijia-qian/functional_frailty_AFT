# ---- Packages ----
library(truncnorm)
library(splines)

# ---- Inverse-Gamma sampler (shape A, scale B): p(x) ∝ x^{-(A+1)} exp(-B/x) ----
rinv_gamma_shape_scale <- function(n, shape, scale) {
  1 / rgamma(n, shape = shape, rate = scale)
}

# ---- Penalty matrix D (PD if a > 0) ----
penalty_matrix <- function(kp, nS, a) {
  D <- nS
  s <- seq(0, 1, length.out = nS)
  spline_basis <- bs(s, df = kp, intercept = TRUE)  # nS x kp
  
  diff2 <- matrix(
    rep(c(1, -2, 1, rep(0, D - 3)), D - 2)[1:((D - 2) * D)],
    nrow = D - 2, ncol = D, byrow = TRUE
  )
  
  P2  <- t(spline_basis) %*% t(diff2) %*% diff2 %*% spline_basis
  Pen <- a * diag(kp) + (1 - a) * P2
  Pen
}

# ---- Trapezoid integration weights on a grid ----
trapz_weights <- function(s_grid) {
  P <- length(s_grid)
  ds <- diff(s_grid)
  w <- numeric(P)
  w[1] <- ds[1] / 2
  w[P] <- ds[P - 1] / 2
  if (P > 2) w[2:(P - 1)] <- (ds[1:(P - 2)] + ds[2:(P - 1)]) / 2
  w
}

# ---- Sample MVN given PRECISION Q and linear term h ----
# If p(theta|.) ∝ exp(-1/2 theta^T Q theta + theta^T h), then
# theta ~ N(Q^{-1} h, Q^{-1})
rmvn_precision <- function(Q, h) {
  # Q must be symmetric positive definite
  L <- chol(Q)  # upper triangular, Q = t(L) %*% L
  mu <- backsolve(L, forwardsolve(t(L), h))
  z <- rnorm(length(h))
  v <- backsolve(L, z)          # Cov(v)=Q^{-1}
  as.vector(mu + v)
}

# ---- MVN draw given covariance Cholesky ----
rmvn_chol <- function(mu, chol_V) {
  mu + as.vector(t(chol_V) %*% rnorm(length(mu)))
}

# ---- Gibbs sampler ----
gibbs_functional_frailty <- function( 
    
    time, status, cluster_id, Z, X, s_grid,
    K = 12,
    a_pen = 0.001,
    lambda = 1000,
    sigma_gamma2 = 25,
    A = 2, B = 1,
    n_iter = 4000,
    n_burn = 1000,
    n_thin = 1,
    seed = 123,
    verbose = TRUE
) {
  set.seed(seed)
  
  N <- length(time)
  stopifnot(length(status) == N, length(cluster_id) == N)
  stopifnot(nrow(Z) == N, nrow(X) == N)
  P <- length(s_grid)
  stopifnot(ncol(X) == P)
  
  cluster_id <- as.integer(factor(cluster_id))
  J <- max(cluster_id)
  
  y_obs <- log(time)
  c_log <- y_obs
  cens_idx <- which(status == 0L)
  
  # build W = ∫ X(s) phi(s) ds with phi from bs(s_grid, df=K)
  Bmat <- bs(s_grid, df = K, intercept = TRUE)  # P x K
  w_int <- trapz_weights(s_grid)
  Xw <- sweep(X, 2, w_int, `*`)
  W <- Xw %*% Bmat  # N x K
  
  # # assessment indicator
  # is_post1 <- df_wide$mtm == "1_Post1"
  # is_post2 <- df_wide$mtm == "2_Post2"
  # stopifnot(all(is_post1 + is_post2 == 1L))  # each row must be exactly one assessment
  # W1 <- W0 * as.numeric(is_post1)  # multiplies each row by 0/1
  # W2 <- W0 * as.numeric(is_post2)
  # W <- cbind(W1, W2)          # N x (2K)
  # K2 <- 2 * K
  
  # penalty matrix D (K x K)
  Dmat <- penalty_matrix(kp = K, nS = P, a = a_pen)
  # # block-diagonal penalty D_big
  # Dmat <- matrix(0, nrow = K2, ncol = K2)
  # Dmat[1:K, 1:K] <- D0
  # Dmat[(K + 1):K2, (K + 1):K2] <- D0
  
  # cluster indices
  idx_by_cluster <- split(seq_len(N), cluster_id)
  
  # init
  M <- ncol(Z)
  gamma <- rep(0, M)
  b <- rep(0, K)
  u <- rep(0, J)
  sigma2 <- max(var(y_obs), 1e-3)
  tau2 <- 1
  
  y_star <- y_obs
  # ensure y* > log(c_ij) for censored
  if (length(cens_idx) > 0) y_star[cens_idx] <- y_obs[cens_idx] + abs(rnorm(length(cens_idx), 0, 0.1)) 
  
  # storage
  keep_idx <- seq(from = n_burn + 1, to = n_iter, by = n_thin)
  S <- length(keep_idx)
  out <- list(
    gamma = matrix(NA_real_, S, M),
    b     = matrix(NA_real_, S, K),
    u     = matrix(NA_real_, S, J),
    sigma2 = numeric(S),
    tau2   = numeric(S),
    beta_mean = NULL,
    beta_q025 = NULL,
    beta_q975 = NULL,
    meta = list(N = N, J = J, M = M, K = K, P = P, lambda = lambda, a_pen = a_pen,
                A = A, B = B, sigma_gamma2 = sigma_gamma2,
                n_iter = n_iter, n_burn = n_burn, n_thin = n_thin)
  )
  colnames(out$gamma) <- colnames(Z)
  
  # # precompute crossprods
  # ZtZ <- crossprod(Z)
  # WtW <- crossprod(W)
  
  # Precompute X_gb = [Z W] and its crossproduct
  Xgb <- cbind(Z, W)          # N x (M+K)
  XtX_gb <- crossprod(Xgb)    # (M+K) x (M+K)
  
  # Prior precision for (gamma,b): blockdiag(1/sigma_gamma2 I, lambda D)
  Q0_gb <- matrix(0, nrow = M + K, ncol = M + K)
  Q0_gb[1:M, 1:M] <- diag(1 / sigma_gamma2, M)
  Q0_gb[(M + 1):(M + K), (M + 1):(M + K)] <- lambda * Dmat
  
  save_counter <- 0L
  for (it in seq_len(n_iter)) {
    
    # 0) impute y* for censored
    mu_all <- as.vector(Z %*% gamma + W %*% b + u[cluster_id])
    if (length(cens_idx) > 0) {
      y_star[cens_idx] <- truncnorm::rtruncnorm(
        n = length(cens_idx),
        a = c_log[cens_idx],
        b = Inf,
        mean = mu_all[cens_idx],
        sd = sqrt(sigma2)
      )
    }
    
    # 1) Block update for theta_gb = (gamma, b) | u, sigma2, y_star
    # Model: y_tilde = y_star - u[cluster_id] = Xgb * theta_gb + eps
    y_tilde <- y_star - u[cluster_id]
    
    Q_gb <- (XtX_gb / sigma2) + Q0_gb
    h_gb <- as.vector(crossprod(Xgb, y_tilde) / sigma2)
    
    theta_gb <- rmvn_precision(Q_gb, h_gb)
    gamma <- theta_gb[1:M]
    b     <- theta_gb[(M + 1):(M + K)]
    
    # # 1) gamma | rest
    # eg <- y_star - as.vector(W %*% b) - u[cluster_id]
    # Qg <- (ZtZ / sigma2) + diag(1 / sigma_gamma2, M)
    # Vg <- solve(Qg)
    # mug <- Vg %*% (crossprod(Z, eg) / sigma2)
    # gamma <- rmvn_chol(as.vector(mug), chol(Vg))
    # 
    # # 2) b | rest
    # eb <- y_star - as.vector(Z %*% gamma) - u[cluster_id]
    # Qb <- (WtW / sigma2) + lambda * Dmat
    # Vb <- solve(Qb)
    # mub <- Vb %*% (crossprod(W, eb) / sigma2)
    # b <- rmvn_chol(as.vector(mub), chol(Vb))
    
    # 3) u_j | rest (independent)
    eu <- y_star - as.vector(Z %*% gamma) - as.vector(W %*% b)
    for (j in seq_len(J)) {
      idx <- idx_by_cluster[[j]]
      Vuj <- 1 / (length(idx) / sigma2 + 1 / tau2)
      muj <- Vuj * (sum(eu[idx]) / sigma2)
      u[j] <- rnorm(1, muj, sqrt(Vuj))
    }
    
    # 4) sigma2 | rest
    r <- y_star - as.vector(Z %*% gamma) - as.vector(W %*% b) - u[cluster_id]
    SSE <- sum(r^2)
    sigma2 <- rinv_gamma_shape_scale(1, shape = A + N / 2, scale = B + 0.5 * SSE)
    
    # 5) tau2 | rest
    tau2 <- rinv_gamma_shape_scale(1, shape = A + J / 2, scale = B + 0.5 * sum(u^2))
    
    # save
    if (it %in% keep_idx) {
      save_counter <- save_counter + 1L
      out$gamma[save_counter, ] <- gamma
      out$b[save_counter, ] <- b
      out$u[save_counter, ] <- u
      out$sigma2[save_counter] <- sigma2
      out$tau2[save_counter] <- tau2
    }
    
    if (verbose && (it %% 1000 == 0)) {
      cat(sprintf("iter %d/%d | sigma2=%.4f tau2=%.4f\n", it, n_iter, sigma2, tau2))
    }
  }
  
  # Posterior summaries of beta(s) on s_grid
  beta_draws <- out$b %*% t(Bmat)  # S x P
  out$beta_mean <- colMeans(beta_draws)
  out$beta_q025 <- apply(beta_draws, 2, quantile, 0.025)
  out$beta_q975 <- apply(beta_draws, 2, quantile, 0.975)
  
  out
}


