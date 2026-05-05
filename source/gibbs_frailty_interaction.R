gibbs_frailty_interaction <- function( 
    time, status, cluster_id, Z, X, age, s_grid,
    K_s = 10,         # Basis dimension for Time (MIMS)
    K_age = 5,        # Basis dimension for Age
    a_pen = 1e-3,     # Small ridge penalty for tensor null-space stability
    bs_s = "cc",    # "cc" for real NHANES data, "cr" for simulation
    bs_age = "cr",
    # tuning / priors
    lambda_init = 5000,
    A_lambda = 1, B_lambda = 0.001,
    var_gamma = 100,
    A_tau2 = 3, B_tau2 = 2,
    A_sigma2 = 3, B_sigma2 = 2,
    # MCMC
    n_iter = 4000,
    n_burn = 1000,
    n_thin = 1,
    verbose = TRUE
) {
  
  require(mgcv)
  require(truncnorm)
  
  N <- length(time)
  stopifnot(length(status) == N, length(cluster_id) == N, length(age) == N)
  stopifnot(nrow(Z) == N, nrow(X) == N)
  P <- length(s_grid)
  stopifnot(ncol(X) == P)
  
  # ---- 1. Prepare Cluster IDs ----
  cluster_factor <- factor(cluster_id)
  cluster_levels <- levels(cluster_factor)
  cluster_id <- as.integer(cluster_factor)
  J <- max(cluster_id)
  
  y_obs <- log(time)
  c_log <- y_obs
  cens_idx <- which(status == 0L)
  
  # ---- 2. Construct Tensor Product Basis & Penalties ----
  # Dummy data spanning the ranges to set knots correctly
  dummy_data <- data.frame(
    s_grid = rep(s_grid, N),
    age    = rep(age, each = P)
  )
  
  # cc = cyclic cubic (for 24h wrap-around), cr = cubic regression (for age)
  sm_tensor <- smoothCon(
    te(s_grid, age, bs = c(bs_s, bs_age), k = c(K_s, K_age)),
    data = dummy_data, absorb.cons = FALSE
  )[[1]]
  
  S1 <- sm_tensor$S[[1]] # Penalty for the 'cc' marginal (time)
  S2 <- sm_tensor$S[[2]] # Penalty for the 'cr' marginal (age)
  
  # rank_S1 <- qr(S1)$rank
  # rank_S2 <- qr(S2)$rank
  rank_S1 <- sm_tensor$rank[1]
  rank_S2 <- sm_tensor$rank[2]
  K_total <- ncol(sm_tensor$X)
  
  # ---- 3. Compute the Integrated Design Matrix W ----
  if (verbose) cat("Pre-computing integrated design matrix W...\n")
  W <- matrix(0, nrow = N, ncol = K_total)
  w_dt <- trapz_weights(s_grid) 
  
  for (i in 1:N) {
    df_i <- data.frame(s_grid = s_grid, age = rep(age[i], P))
    Phi_i <- PredictMat(sm_tensor, df_i) # Size: P x K_total
    
    # Integrate: sum_t( X_i(t) * Phi_i(t, age_i) * w_dt )
    W[i, ] <- colSums(X[i, ] * Phi_i * w_dt)
  }
  
  # ---- 4. MCMC Initialization ----
  idx_by_cluster <- split(seq_len(N), cluster_id)
  
  M <- ncol(Z)
  gamma  <- rep(0, M)
  b      <- rep(0, K_total)
  u      <- rep(0, J)
  sigma2 <- 1
  tau2   <- 1
  lambda1 <- lambda_init
  lambda2 <- lambda_init
  
  y_star <- y_obs
  if (length(cens_idx) > 0) {
    y_star[cens_idx] <- y_obs[cens_idx] + abs(rnorm(length(cens_idx), 0, 0.1)) 
  }
  
  # Storage setup
  keep_idx <- seq(from = n_burn + 1, to = n_iter, by = n_thin)
  S <- length(keep_idx)
  out <- list(
    gamma  = matrix(NA_real_, S, M),
    b      = matrix(NA_real_, S, K_total),
    u      = matrix(NA_real_, S, J),
    sigma2 = numeric(S),
    tau2   = numeric(S),
    lambda1 = numeric(S),
    lambda2 = numeric(S),
    meta = list(
      N = N, J = J, M = M, K_total = K_total, P = P, a_pen = a_pen, 
      s_grid = s_grid, cluster_levels = cluster_levels, 
      rank_S1 = rank_S1, rank_S2 = rank_S2,
      n_iter = n_iter, n_burn = n_burn, n_thin = n_thin
    )
  )
  if (!is.null(colnames(Z))) colnames(out$gamma) <- colnames(Z)
  
  Xgb <- cbind(Z, W)          
  XtX_gb <- crossprod(Xgb)
  
  idx_b <- (M + 1):(M + K_total)
  diag_ridge <- diag(a_pen, K_total)
  
  Q0_gb <- matrix(0, nrow = M + K_total, ncol = M + K_total)
  diag(Q0_gb)[1:M] <- 1 / var_gamma
  
  # ---- 5. Gibbs Sampler Loop ----
  save_counter <- 0L
  if (verbose) cat("Starting MCMC...\n")
  
  for (it in seq_len(n_iter)) {
    
    # 1) Impute y* for censored subjects
    mu_all <- as.vector(Xgb %*% c(gamma, b) + u[cluster_id])
    if (length(cens_idx) > 0) {
      y_star[cens_idx] <- rtruncnorm(
        n = length(cens_idx),
        a = c_log[cens_idx],
        b = Inf,
        mean = mu_all[cens_idx],
        sd = sqrt(sigma2)
      )
    }
    
    # 2) Block update for theta_gb = (gamma, b) | u, sigma2, y_star
    y_tilde <- y_star - u[cluster_id]
    
    # Combine the two penalties for the block
    Q0_gb[idx_b, idx_b] <- lambda1 * S1 + lambda2 * S2 + diag_ridge
    
    Q_gb <- (XtX_gb / sigma2) + Q0_gb
    h_gb <- as.vector(crossprod(Xgb, y_tilde) / sigma2)
    
    theta_gb <- rmvn_precision(Q_gb, h_gb)
    gamma <- theta_gb[1:M]
    b     <- theta_gb[idx_b]
    
    # 3) Update frailties u_j | rest
    eu <- y_star - as.vector(Xgb %*% c(gamma, b))
    for (j in seq_len(J)) {
      idx <- idx_by_cluster[[j]]
      Vuj <- 1 / (length(idx) / sigma2 + 1 / tau2)
      muj <- Vuj * (sum(eu[idx]) / sigma2)
      u[j] <- rnorm(1, muj, sqrt(Vuj))
    }
    
    # 4) Variances: sigma2 (residuals) and tau2 (frailty) | rest
    r <- y_star - as.vector(Xgb %*% c(gamma, b)) - u[cluster_id]
    SSE <- sum(r^2)
    sigma2 <- rinv_gamma_shape_scale(1, shape = A_sigma2 + N / 2, scale = B_sigma2 + 0.5 * SSE)
    
    tau2 <- rinv_gamma_shape_scale(1, shape = A_tau2 + J / 2, scale = B_tau2 + 0.5 * sum(u^2))
    
    # 5) Smoothing parameters: lambda1 (time) and lambda2 (age) | b
    quad_b1 <- as.numeric(t(b) %*% S1 %*% b)
    lambda1 <- rgamma(1, shape = A_lambda + rank_S1 / 2, rate = B_lambda + 0.5 * quad_b1)
    
    quad_b2 <- as.numeric(t(b) %*% S2 %*% b)
    lambda2 <- rgamma(1, shape = A_lambda + rank_S2 / 2, rate = B_lambda + 0.5 * quad_b2)
    
    # ---- Save & Print ----
    if (it %in% keep_idx) {
      save_counter <- save_counter + 1L
      out$gamma[save_counter, ] <- gamma
      out$b[save_counter, ]     <- b
      out$u[save_counter, ]     <- u
      out$sigma2[save_counter]  <- sigma2
      out$tau2[save_counter]    <- tau2
      out$lambda1[save_counter] <- lambda1
      out$lambda2[save_counter] <- lambda2
    }
    
    if (verbose && it %% 1000 == 0) {
      cat(sprintf(
        "iter %d/%d | sigma2=%.4f tau2=%.4f lam1=%.2f lam2=%.2f\n",
        it, n_iter, sigma2, tau2, lambda1, lambda2
      ))
    }
  }
  
  # Note: The calculation of the 3D surface (beta_mean, etc.) requires mapping 
  # the posterior 'b' draws back through a 2D prediction matrix. 
  # This is usually best done outside the main Gibbs loop to save memory.
  
  # Append the smoothCon object to the meta list so you can use it for predictions later
  out$meta$sm_tensor <- sm_tensor 
  
  return(out)
}