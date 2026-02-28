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

# ---- CMA band from beta_draws (Q x P) ----
cma_band <- function(beta_draws, level = 0.95, eps = 1e-12) {
  stopifnot(is.matrix(beta_draws))
  Q <- nrow(beta_draws); P <- ncol(beta_draws)
  
  beta_hat <- colMeans(beta_draws)
  S_hat <- apply(beta_draws, 2, sd)
  S_use <- pmax(S_hat, eps)
  
  # d^q = max_t |beta^q(t)-beta_hat(t)| / S_hat(t)
  d <- apply(
    abs(beta_draws - matrix(beta_hat, Q, P, byrow = TRUE)) /
      matrix(S_use, Q, P, byrow = TRUE),
    1, max
  )
  
  # critical value
  qd <- as.numeric(quantile(d, probs = level, names = FALSE))
  
  list(
    mean = beta_hat,
    sd = S_hat,
    d = d,
    qd = qd,
    lower = beta_hat - qd * S_hat,
    upper = beta_hat + qd * S_hat
  )
}

# ---- Gibbs sampler ----
gibbs_functional_frailty <- function( 
    time, status, cluster_id, Z, X, s_grid,
    K = 12,
    a_pen = 0.001,
    
    # lambda prior + init (Gamma shape-rate)
    lambda_init = 1000,
    A_lambda = 1,
    B_lambda = 0.001,
    lambda = 1000,
    
    var_gamma = 100,
    A_tau2 = 3, B_tau2 = 2,
    A_sigma2 = 3, B_sigma2 = 2,
    n_iter = 4000,
    n_burn = 1000,
    n_thin = 1,
    verbose = TRUE
) {
  
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
  
  # penalty matrix D (K x K), PD if a_pen > 0
  Dmat <- penalty_matrix(kp = K, nS = P, a = a_pen)
  
  # cluster indices
  idx_by_cluster <- split(seq_len(N), cluster_id)
  
  # init
  M <- ncol(Z)
  gamma  <- rep(0, M)
  b      <- rep(0, K)
  u      <- rep(0, J)
  sigma2 <- 1
  tau2   <- 1
  lambda <- lambda_init
  
  y_star <- y_obs
  if (length(cens_idx) > 0) {
    y_star[cens_idx] <- y_obs[cens_idx] + abs(rnorm(length(cens_idx), 0, 0.1)) 
  }
  
  # storage
  keep_idx <- seq(from = n_burn + 1, to = n_iter, by = n_thin)
  S <- length(keep_idx)
  out <- list(
    gamma  = matrix(NA_real_, S, M),
    b      = matrix(NA_real_, S, K),
    u      = matrix(NA_real_, S, J),
    sigma2 = numeric(S),
    tau2   = numeric(S),
    lambda = numeric(S),
    beta_mean = NULL,
    beta_q025 = NULL,
    beta_q975 = NULL,
    beta_cma_lower = NULL,
    beta_cma_upper = NULL,
    beta_cma_qd    = NULL,
    beta_cma_sd    = NULL,
    meta = list(
      N = N, J = J, M = M, K = K, P = P, a_pen = a_pen, 
      lambda_init = lambda_init,
      A_lambda = A_lambda, B_lambda = B_lambda,
      A_sigma2 = A_sigma2, A_tau2 = A_tau2, B_sigma2 = B_sigma2, B_tau2 = B_tau2,
      var_gamma = var_gamma,
      n_iter = n_iter, n_burn = n_burn, n_thin = n_thin
    )
  )
  colnames(out$gamma) <- colnames(Z)
  
  # Precompute X_gb = [Z W] and its crossproduct
  Xgb <- cbind(Z, W)          
  XtX_gb <- crossprod(Xgb)
  
  # ---- IMPROVEMENT: preallocate Q0_gb once, only update b-block ----
  idx_b <- (M + 1):(M + K)
  
  Q0_gb <- matrix(0, nrow = M + K, ncol = M + K)
  Q0_gb[1:M, 1:M] <- diag(1 / var_gamma, M)
  # Ensure off-diagonal blocks are zero once (they already are):
  # Q0_gb[1:M, idx_b] <- 0
  # Q0_gb[idx_b, 1:M] <- 0
  # Initialize b-block:
  Q0_gb[idx_b, idx_b] <- lambda * Dmat
  
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
    y_tilde <- y_star - u[cluster_id]
    
    # update only the KxK block with current lambda
    Q0_gb[idx_b, idx_b] <- lambda * Dmat
    
    Q_gb <- (XtX_gb / sigma2) + Q0_gb
    h_gb <- as.vector(crossprod(Xgb, y_tilde) / sigma2)
    
    theta_gb <- rmvn_precision(Q_gb, h_gb)
    gamma <- theta_gb[1:M]
    b     <- theta_gb[idx_b]
    
    # 3) u_j | rest
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
    sigma2 <- rinv_gamma_shape_scale(
      1,
      shape = A_sigma2 + N / 2,
      scale = B_sigma2 + 0.5 * SSE
    )
    
    # 5) tau2 | rest
    tau2 <- rinv_gamma_shape_scale(
      1,
      shape = A_tau2 + J / 2,
      scale = B_tau2 + 0.5 * sum(u^2)
    )
    
    # 6) lambda | b  (Gamma full conditional; shape-rate)
    quad_b <- as.numeric(t(b) %*% Dmat %*% b)
    lambda <- rgamma(
      n = 1,
      shape = A_lambda + K / 2,
      rate  = B_lambda + 0.5 * quad_b
    )
    # lambda <- 10000
    
    # save
    if (it %in% keep_idx) {
      save_counter <- save_counter + 1L
      out$gamma[save_counter, ] <- gamma
      out$b[save_counter, ]     <- b
      out$u[save_counter, ]     <- u
      out$sigma2[save_counter]  <- sigma2
      out$tau2[save_counter]    <- tau2
      out$lambda[save_counter]  <- lambda
    }
    
    if (verbose && it %% 1000 == 0) {
      cat(sprintf(
        "iter %d/%d | sigma2=%.4f tau2=%.4f lambda=%.2f\n",
        it, n_iter, sigma2, tau2, lambda
      ))
    }
  }
  
  # Posterior summaries of beta(s) on s_grid
  beta_draws <- out$b %*% t(Bmat)  # S x P
  out$beta_mean <- colMeans(beta_draws)
  out$beta_q025 <- apply(beta_draws, 2, quantile, 0.025)
  out$beta_q975 <- apply(beta_draws, 2, quantile, 0.975)
  
  # CMA band
  cma <- cma_band(beta_draws, level = 0.95, eps = 1e-12)
  out$beta_cma_lower <- cma$lower
  out$beta_cma_upper <- cma$upper
  out$beta_cma_qd    <- cma$qd
  out$beta_cma_sd    <- cma$sd
  
  out
}
