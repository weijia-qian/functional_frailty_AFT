simulate_AFT = function(data = dat_func, 
                        family = c("lognormal", "loglogistic"),
                        n_cluster = 100,
                        n_subject = 5,
                        # npc = 5,
                        tmax = 1,
                        nS = 401,
                        k0 = 20,
                        alpha = 0.7,
                        beta_type = c("monotone", "peak1", "peak2", "wavy"),
                        gamma = c(0.5, 0.3, -0.2), # first is intercept
                        sigma = 0.2,
                        tau = 0.5,
                        censor_rate = 0.25
                        # seed = NULL
                        ){
  
  family <- match.arg(family)
  beta_type <- match.arg(beta_type)
  
  # ---- choose beta basis coefficients ----
  if (beta_type == "monotone") {
    bs_coef = c(0.1, -0.15, -0.35, -0.5, -0.5, -0.5)
  } else if (beta_type == "peak1") {
    bs_coef = c(0.0, -0.1, 0.4, 1.2, 0.5, 0.1, 0.0)
  } else if (beta_type == "peak2") {
    bs_coef = c(0.0, -0.1, 0.4, 0.9, 0.3, -0.15, -0.3, 0.2, 0.8, 0.3, 0.0)
  } else if (beta_type == "wavy") {
    bs_coef = c(0.0, 1.7, -2.1, 2.5, -0.85, 0.6, -0.4, 0.2, 0.0)
    # bs_coef = c(0.0, 1.7, -2.1, 2.5, -2.3, 2.1, -1.9, 1.5, -0.85, 0.6, -0.4, 0.2, 0.0)
    # bs_coef = c(0.0, 1.7, -2.1, 2.5, -2.3, 2.1, -1.9, 1.5, -1.25, 1.1, -0.85, 0.6, -0.4, 0.2, 0.0)
  }
  
  # # ---- seed ----
  # if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  # set.seed(seed)
  
  # ---- build cluster/subject indexing ----
  N = n_cluster * n_subject
  cluster_id <- rep(seq_len(n_cluster), each = n_subject)
  
  # ---- simulate two scalar covariates (z1, z2) ----
  z1 <- rnorm(N, mean = 0, sd = 1)
  z2 <- rbinom(N, size = 1, prob = 0.5)
  
  # ---- simulation grid ----
  s_grid <- seq(0, tmax, length.out = nS)
  
  # ---- Construct Fourier basis ----
  Phi <- matrix(0, nS, k0)
  
  Phi[,1] <- 1
  
  for (k in 1:floor((k0-1)/2)) {
    Phi[,2*k]   <- sqrt(2) * cos(2*pi*k*s_grid)
    Phi[,2*k+1] <- sqrt(2) * sin(2*pi*k*s_grid)
  }
  
  # ---- Eigenvalue decay ----
  lambda <- (1:k0)^(-alpha)
  
  # ---- Random scores ----
  Xi <- matrix(rnorm(N*k0), N, k0)
  
  # ---- Generate curves ----
  sim_curves <- Xi %*% diag(sqrt(lambda)) %*% t(Phi)

  # # ---- simulate X(s): 0.5 sin(s), 0.5 cos(s) ----
  # # Each subject draws one of the two shapes
  # is_sin <- rbinom(N, size = 1, prob = 0.5)
  # 
  # A_i   <- rnorm(N, mean = 1, sd = 0.2)        # amplitude
  # phi_i <- rnorm(N, mean = 0, sd = 0.15)       # phase shift
  # 
  # base_curves <- matrix(NA_real_, N, nS)
  # for (i in 1:N) {
  #   if (is_sin[i]) base_curves[i,] <- sin(s_grid + phi_i[i])
  #   else           base_curves[i,] <- cos(s_grid + phi_i[i])
  # }
  # base_curves <- base_curves * A_i
  # 
  # # random error curves
  # Qe <- 15
  # Be <- splines::bs(s_grid, df = Qe, intercept = TRUE)   # nS x Qe
  # alpha <- matrix(rnorm(N * Qe, 0, 0.05), N, Qe)         # controls noise size
  # err_curves <- alpha %*% t(Be)                          # N x nS smooth
  # sim_curves <- base_curves + err_curves

  # # ---- FPCA on real data ----
  # df_wide <- data |>
  #   dplyr::select(-dplyr::any_of("seconds")) |>
  #   tidyr::pivot_wider(
  #     names_from = frame,
  #     values_from = percent_change,
  #     names_prefix = "pct_chg_"
  #   )
  # 
  # s_orig <- unique(data$seconds)   # original grid
  # s_grid <- seq(0, tmax, length.out = nS)  # simulation grid
  # 
  # X_cols <- grep("^pct_chg_", names(df_wide), value = TRUE)
  # mat_real <- as.matrix(df_wide[, X_cols, drop = FALSE])
  # 
  # # align / interpolate real curves onto s_grid
  # if (nS == length(s_orig)) {
  #   mat_fpca <- mat_real
  # } else {
  #   # interpolate each row onto nS points using index scale
  #   x_in <- seq(0, tmax, length.out = ncol(mat_real))
  #   x_out <- s_grid
  #   mat_fpca <- t(apply(mat_real, 1, function(row) {
  #     stats::approx(x = x_in, y = row, xout = x_out, rule = 2)$y
  #   }))
  # }
  # 
  # # apply FPCA
  # fpca.results <- refund::fpca.face(mat_fpca, argvals = s_grid, pve = 0.999)
  # 
  # # re-scaled to appear on the scale of the original functions
  # Phi <- sqrt(nS) * fpca.results$efunctions
  # eigenvalues <- fpca.results$evalues / nS
  # 
  # # ---- simulate FPCA scores and curves for N observations ----
  # sim_scores <- mvrnorm(n = N, mu = rep(0, npc), Sigma = diag(eigenvalues[1:npc]))
  # # mu <- matrix(rep(fpca.results$mu, N), nrow = N, byrow = TRUE)
  # mu <- matrix(fpca.results$mu, nrow = N, ncol = nS, byrow = TRUE)
  # sim_curves <- mu + sim_scores %*% t(Phi[, 1:npc])
  # sim_curves[, 1] <- 0 # assign the initial value to zero
  
  # ---- define spline basis for beta(s) on s_grid ----
  k <- length(bs_coef)
  B <- splines::bs(s_grid, df = k, intercept = TRUE) # ns x k
  beta <- B %*% bs_coef
  
  # ---- numerical integration: ∫ X(s) beta(s) ds ----
  wts <- diag(trapz_weights(s_grid)) 
  num_int <- as.vector((sim_curves %*% wts) %*% beta)
  
  Xw <- sim_curves %*% wts
  W <- Xw %*% B  # N x K
  
  # ---- simulate error term z ----
  if (family == "loglogistic"){
    z <- stats::rlogis(N)
  } else if (family == "lognormal"){
    z <- rnorm(N)
    # # use Box-Muller transform
    # u1 <- runif(N)
    # u2 <- runif(N)
    # z <- sqrt(-2 * log(u1)) * cos(2 * pi * u2)
  } else {
    stop('Invalid family.')
  }
  
  # ---- shared frailty u_j by subject ----
  u <- stats::rnorm(n_cluster, mean = 0, sd = tau)
  frailty <- u[cluster_id]
  
  # ---- generate survival time ----
  lp <- gamma[1] + gamma[2] * z1 + gamma[3] * z2 + num_int
  T_true <- exp(lp + frailty + sigma * z)
  
  # ---- censoring ----
  if (censor_rate > 0){
    ub <- choose_ub_unif(T_true, censor_rate)
    C <- stats::runif(N, 0, ub)
    Y <- pmin(T_true, C)
  } else {
    Y <- T_true
  }
  delta <- as.integer(Y == T_true)
  
  # ---- save simulated data ----
  sim_data = data.frame(subject_id = seq(1:N), 
                        cluster_id = cluster_id,
                        Y = Y,
                        delta = delta,
                        X = I(sim_curves),
                        W = I(W),
                        Z1 = z1,
                        Z2 = z2,
                        T_true = T_true,
                        # C = C,
                        lp = lp)
  
  # ---- true coefficient functions ----
  df_coef = data.frame(time = s_grid,
                       beta = beta,
                       gamma = I(matrix(gamma, ncol = length(gamma), nrow = nS, byrow = TRUE)),
                       sigma = rep(sigma, nS),
                       tau = rep(tau, nS)
                       )
  
  out <- list(data = sim_data, coefficients = df_coef, family = family, beta_type = beta_type, bs_coef = bs_coef, u = u)
  
  out
}

trapz_weights <- function(x) {
  n <- length(x)
  dx <- diff(x)
  w <- numeric(n)
  w[1] <- dx[1] / 2
  w[n] <- dx[n - 1] / 2
  if (n > 2) w[2:(n - 1)] <- (dx[1:(n - 2)] + dx[2:(n - 1)]) / 2
  w
}

choose_ub_unif <- function(T_true, censor_rate, tol = 1e-8) {
  stopifnot(censor_rate > 0, censor_rate < 1)
  
  g <- function(ub) mean(pmin(T_true, ub) / ub) - censor_rate
  
  # Lower/upper bracket for uniroot
  lo <- max(tol, min(T_true[T_true > 0], na.rm = TRUE) * 1e-6)
  hi <- max(T_true, na.rm = TRUE) * 2 + 1
  
  # Ensure bracket has opposite signs: g(lo) ~ 1 - censor_rate > 0, g(hi) < 0
  if (g(lo) < 0) lo <- tol
  while (g(hi) > 0) hi <- hi * 2
  
  uniroot(g, lower = lo, upper = hi)$root
}


