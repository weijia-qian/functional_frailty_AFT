simulate_AFT = function(data = dat_func, 
                        family = c("lognormal", "loglogistic"),
                        n_cluster = 100,
                        n_subject = 5,
                        npc = 5,
                        tmax = 1,
                        nS = 401,
                        beta_type = c("simple", "complex"),
                        beta0 = 0.5,
                        sigma = 0.2,
                        tau = 0.5,
                        gamma = c(0.3, -0.2),
                        seed = NULL
                        ){
  
  family <- match.arg(family)
  beta_type <- match.arg(beta_type)
  
  # ---- choose beta basis coefficients + censoring upper bound ----
  if (beta_type == "simple") {
    bs_coef = c(0, -1, -0.5, 0.25, 0.25, 0.25)
    ub = 250
  } else if (beta_type == "complex") {
    bs_coef = c(0, -0.6, -1.2, 0.6, -0.5, 1, 0.5, 0)
    ub = 35
  }
  
  # ---- seed ----
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  set.seed(seed)
  
  # ---- build cluster/subject indexing ----
  N = n_cluster * n_subject
  cluster_id <- rep(seq_len(n_cluster), each = n_subject)
  
  # ---- simulate two scalar covariates (z1, z2) ----
  # z1 ~ N(0,1), z2 ~ Bernoulli(0.5) by default
  z1 <- rnorm(N, mean = 0, sd = 1)
  z2 <- rbinom(N, size = 1, prob = 0.5)
  
  # ---- FPCA on real data ----
  df_wide <- data |>
    dplyr::select(-dplyr::any_of("seconds")) |>
    tidyr::pivot_wider(
      names_from = frame,
      values_from = percent_change,
      names_prefix = "pct_chg_"
    )
  
  s_orig <- unique(data$seconds)   # original grid
  s_grid <- seq(0, tmax, length.out = nS)  # simulation grid
  
  X_cols <- grep("^pct_chg_", names(df_wide), value = TRUE)
  mat_real <- as.matrix(df_wide[, X_cols, drop = FALSE])
  
  # align / interpolate real curves onto s_grid
  if (nS == length(s_orig)) {
    mat_fpca <- mat_real
  } else {
    # interpolate each row onto nS points using index scale
    x_in <- seq(0, tmax, length.out = ncol(mat_real))
    x_out <- s_grid
    mat_fpca <- t(apply(mat_real, 1, function(row) {
      stats::approx(x = x_in, y = row, xout = x_out, rule = 2)$y
    }))
  }
  
  # apply FPCA
  fpca.results <- refund::fpca.face(mat_fpca, argvals = s_grid, pve = 0.99)
  
  # re-scaled to appear on the scale of the original functions 
  Phi <- sqrt(nS) * fpca.results$efunctions
  eigenvalues <- fpca.results$evalues / nS
  
  # ---- simulate FPCA scores and curves for N observations ----
  sim_scores <- mvrnorm(n = N, mu = rep(0, npc), Sigma = diag(eigenvalues[1:npc]))
  # mu <- matrix(rep(fpca.results$mu, N), nrow = N, byrow = TRUE)
  mu <- matrix(fpca.results$mu, nrow = N, ncol = nS, byrow = TRUE)
  sim_curves <- mu + sim_scores %*% t(Phi[, 1:npc])
  sim_curves[, 1] <- 0 # assign the initial value to zero
  
  # ---- define spline basis for beta(s) on s_grid ----
  k <- length(bs_coef)
  B <- splines::bs(s_grid, df = k, intercept = TRUE) # ns x k
  beta1 <- B %*% bs_coef
  
  # ---- numerical integration: ∫ X(s) beta(s) ds ----
  trapz_weights <- function(x) {
    n <- length(x)
    dx <- diff(x)
    w <- numeric(n)
    w[1] <- dx[1] / 2
    w[n] <- dx[n - 1] / 2
    if (n > 2) w[2:(n - 1)] <- (dx[1:(n - 2)] + dx[2:(n - 1)]) / 2
    w
  }
  wts <- trapz_weights(s_grid)
  num_int <- as.vector((sim_curves * wts) %*% beta1)
  
  # ---- simulate error term z ----
  if (family == "loglogistic"){
    z <- stats::rlogis(N)
  } else if (family == "lognormal"){
    #z <- rnorm(n)
    # use Box-Muller transform
    u1 <- runif(N)
    u2 <- runif(N)
    z <- sqrt(-2 * log(u1)) * cos(2 * pi * u2)
  } else {
    stop('Invalid family.')
  }
  
  # ---- shared frailty u_j by subject ----
  u <- stats::rnorm(n_cluster, mean = 0, sd = tau)
  frailty <- u[cluster_id]
  
  # ---- generate survival time ----
  lp <- beta0 + gamma[1] * z1 + gamma[2] * z2 + num_int
  T_true <- exp(lp + frailty + sigma * z)
  
  # ---- censoring ----
  C <- stats::runif(N, 0, ub)
  Y <- pmin(T_true, C)
  delta <- as.integer(T_true <= C)
  
  # ---- save simulated data ----
  sim_data = data.frame(subject_id = seq(1:N), 
                        cluster_id = cluster_id,
                        Y = Y,
                        delta = delta,
                        X = I(sim_curves),
                        Z1 = z1,
                        Z2 = z2,
                        T_true = T_true,
                        C = C,
                        lp = lp)
  
  # ---- true coefficient functions ----
  df_coef = data.frame(time = s_grid,
                       beta0 = rep(beta0, nS),
                       beta1 = beta1,
                       sigma = rep(sigma, nS),
                       tau = rep(tau, nS),
                       gamma1 = rep(gamma[1], nS),
                       gamma2 = rep(gamma[2], nS))
  
  out <- list(data = sim_data, coefficients = df_coef, family = family, beta_type = beta_type, bs_coef = bs_coef, ub = ub, seed = seed)
  
  out
}
