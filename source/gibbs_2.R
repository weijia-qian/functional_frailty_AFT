# ---- Packages ----
library(truncnorm)
library(splines)
library(MASS)
library(MCMCpack)

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

# ---- CMA band from beta_draws (Q x P) ----
cma_band <- function(beta_draws, level = 0.95, eps = 1e-12) {
  stopifnot(is.matrix(beta_draws))
  Q <- nrow(beta_draws)
  P <- ncol(beta_draws)
  
  beta_hat <- colMeans(beta_draws)
  S_hat <- apply(beta_draws, 2, sd)
  S_use <- pmax(S_hat, eps)
  
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
    
    lambda_init = 1000,
    A_lambda = 1, B_lambda = 0.001,
    
    var_gamma = 100,
    A_sigma2 = 2, B_sigma2 = 1,
    A_tau2   = 2, B_tau2   = 1,
    
    n_iter = 4000,
    n_burn = 1000,
    n_thin = 1,
    
    cma_level = 0.95,
    verbose = TRUE
) {
  
  ## ---- dimensions & data ----
  N <- length(time)
  stopifnot(length(status) == N, length(cluster_id) == N)
  stopifnot(nrow(Z) == N, nrow(X) == N)
  
  cluster_id <- as.integer(factor(cluster_id))
  J <- max(cluster_id)
  
  y_obs <- log(time)
  cens_idx <- which(status == 0L)
  
  ## ---- functional design ----
  Bmat <- bs(s_grid, df = K, intercept = TRUE)
  w_int <- trapz_weights(s_grid)
  Xw <- sweep(X, 2, w_int, `*`)
  W <- Xw %*% Bmat
  
  ## ---- penalty ----
  Dmat <- penalty_matrix(kp = K, nS = length(s_grid), a = a_pen)
  
  ## ---- cluster indices ----
  idx_by_cluster <- split(seq_len(N), cluster_id)
  n_j <- lengths(idx_by_cluster)
  
  ## ---- initial values ----
  M <- ncol(Z)
  gamma <- rep(0, M)
  b <- rep(0, K)
  u <- rep(0, J)
  # eta    <- rnorm(J)
  tau2   <- B_tau2 / (A_tau2 - 1) # # mean of IG (if A>1)
  # u      <- sqrt(tau2) * eta
  sigma2 <- B_sigma2 / (A_sigma2 - 1)
  lambda <- lambda_init
  
  y_star <- y_obs
  if (length(cens_idx) > 0) {
    y_star[cens_idx] <- y_obs[cens_idx] + abs(rnorm(length(cens_idx), 0, 0.1))
  }
  
  ## ---- storage ----
  keep_idx <- seq(n_burn + 1, n_iter, by = n_thin)
  S <- length(keep_idx)
  
  out <- list(
    gamma  = matrix(NA_real_, S, M),
    b      = matrix(NA_real_, S, K),
    u      = matrix(NA_real_, S, J),
    sigma2 = numeric(S),
    tau2   = numeric(S),
    lambda = numeric(S)
  )
  colnames(out$gamma) <- colnames(Z)
  
  X_gb <- cbind(Z, W)
  XtX_gb <- crossprod(X_gb)
  
  Q0_base <- matrix(0, M + K, M + K)
  Q0_base[1:M, 1:M] <- diag(1 / var_gamma, M)
  
  save_counter <- 0L
  
  ## ---- Gibbs loop ----
  for (it in seq_len(n_iter)) {
    
    ## ---- 0) impute y* ----
    mu_all <- as.vector(Z %*% gamma + W %*% b + u[cluster_id])
    if (length(cens_idx) > 0) {
      y_star[cens_idx] <- truncnorm::rtruncnorm(
        length(cens_idx),
        a = y_obs[cens_idx],
        b = Inf,
        mean = mu_all[cens_idx],
        sd = sqrt(sigma2)
      )
    }
    
    ## ---- 1) (gamma, b) | rest ----
    y_tilde <- y_star - u[cluster_id]
    
    Q0 <- Q0_base
    Q0[(M + 1):(M + K), (M + 1):(M + K)] <- lambda * Dmat
    
    Q_post  <- XtX_gb / sigma2 + Q0
    V_post  <- solve(Q_post)
    mu_post <- V_post %*% (crossprod(X_gb, y_tilde) / sigma2)
    
    theta <- MASS::mvrnorm(1, as.vector(mu_post), V_post)
    gamma <- theta[1:M]
    b <- theta[(M + 1):(M + K)]
    
    ## ---- 2) u | rest ----
    res_u <- y_star - as.vector(Z %*% gamma) - as.vector(W %*% b)    
    n_j <- lengths(idx_by_cluster)
    sum_res_j <- vapply(idx_by_cluster, function(idx) sum(res_u[idx]), numeric(1))
    
    V_u <- 1 / (n_j / sigma2 + 1 / tau2)
    m_u <- V_u * (sum_res_j / sigma2)
    u <- rnorm(J, mean = m_u, sd = sqrt(V_u))
    
    # ## ---- 2) eta | rest ----
    # res_eta <- y_star - as.vector(Z %*% gamma) - as.vector(W %*% b)
    # 
    # sum_res_j <- vapply(
    #   idx_by_cluster,
    #   function(idx) sum(res_eta[idx]),
    #   numeric(1)
    # )
    # 
    # V_eta <- 1 / (n_j * tau2 / sigma2 + 1)
    # m_eta <- V_eta * (sqrt(tau2) * sum_res_j / sigma2)
    # 
    # eta <- rnorm(J, mean = m_eta, sd = sqrt(V_eta))
    # 
    # ## ---- 3) tau2 | eta ----
    # tau2 <- MCMCpack::rinvgamma(
    #   1,
    #   shape = A_tau2 + J / 2,
    #   scale = B_tau2 + 0.5 * sum(eta^2)
    # )
    # 
    # u <- sqrt(tau2) * eta
    
    ## ---- 3) sigma2 | rest ----
    res <- y_star - as.vector(Z %*% gamma) - as.vector(W %*% b) - u[cluster_id]
    sigma2 <- MCMCpack::rinvgamma(n = 1, shape = A_sigma2 + N / 2, scale = B_sigma2 + 0.5 * sum(res^2))
    
    ## ---- 4) tau2 | rest ----
    tau2 <- MCMCpack::rinvgamma( n = 1, shape = A_tau2 + J / 2, scale = B_tau2 + 0.5 * sum(u^2) )
    
    ## ---- 5) lambda | b ----
    quad_b <- as.numeric(t(b) %*% Dmat %*% b)
    lambda <- rgamma(n = 1, shape = A_lambda + K / 2, rate = B_lambda + 0.5 * quad_b)
    # lambda <- 0
    
    ## ---- save ----
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
  
  ## ---- beta summaries ----
  beta_draws <- out$b %*% t(Bmat)
  out$beta_mean <- colMeans(beta_draws)
  out$beta_q025 <- apply(beta_draws, 2, quantile, 0.025)
  out$beta_q975 <- apply(beta_draws, 2, quantile, 0.975)
  
  cma <- cma_band(beta_draws, level = cma_level)
  out$beta_cma_lower <- cma$lower
  out$beta_cma_upper <- cma$upper
  out$beta_cma_qd    <- cma$qd
  out$beta_cma_sd    <- cma$sd
  
  out
}