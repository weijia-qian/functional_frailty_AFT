# =============================================================================
# simulate_Cox.R
#
# Clustered survival data from a functional Cox (Weibull PH) model with a shared
# Gaussian frailty and an OPTIONAL functional covariate effect:
#
#   h(t_ij | u_j) = h_0(t) * exp( gamma[1] + gamma[2]*Z1_ij + gamma[3]*Z2_ij
#                                 + integral_0^tmax X_ij(s) beta(s) ds + u_j )
#
#   h_0(t) = lambda_0 * rho * t^(rho - 1)     Weibull baseline (rho = 1: exponential)
#   u_j    ~ N(0, tau^2)                      frailty, on the LOG-HAZARD scale
#   T_ij   = ( E_ij / (lambda_0 * exp(lp_ij)) )^(1/rho),   E ~ Exp(1)
#
# This is the beta(s) counterpart of simulate_Cox_interaction.R: the bivariate
# surface beta(s, a) is replaced by a curve, so the functional effect no longer
# varies with age and `age` is not simulated at all.  beta(s) is built from the
# same B-spline coefficient vectors as simulate_AFT.R, so the two describe the
# same coefficient functions -- one on the log-hazard scale, the other on the
# log-time scale.
#
# OPTIONAL FUNCTIONAL EFFECT
# --------------------------
# simulate_X = FALSE omits X(s) and the integral entirely, leaving a Weibull PH
# frailty model in Z1 and Z2 alone.  For a NULL scenario that keeps X -- fit a
# functional model and check that beta_hat(s) ~ 0 -- use beta_scale = 0 instead,
# which simulates the curves but gives them no effect.
#
# CLUSTER DESIGN
# --------------
# N_total is fixed exactly.  All J cluster sizes are drawn from
# Uniform(nj_min, nj_max) with
#   nj_min = floor(mean_nj * (1 - range_factor)),  mean_nj = N_total / n_cluster
#   nj_max = ceil( mean_nj * (1 + range_factor))
# and the discrepancy against N_total is then spread one subject at a time over
# randomly chosen clusters, so every cluster stays close to its nominal
# distribution and no cluster absorbs the whole remainder.  Note this means a
# cluster can end up marginally outside [nj_min, nj_max].
#
# VIOLATING PROPORTIONAL HAZARDS  (psi != 0)
# ------------------------------------------
# With psi = 0 the functional effect is constant in time and the model is both
# PH and (because the baseline is Weibull) AFT.  Setting psi != 0 makes the
# coefficient function drift with log time,
#
#   beta(s, t) = beta(s) * (1 + psi * log t),
#
# so the log hazard ratio between two subjects changes with t and PH fails.
# The hazard becomes
#
#   h(t) = lambda_0 * rho * exp(lp) * t^(rho - 1 + psi*eta),
#          eta = integral X(s) beta(s) ds,
#
# i.e. each subject gets their OWN Weibull shape rho_i = rho + psi*eta_i: the
# functional covariate now modulates the shape, not just the scale.  The
# cumulative hazard still integrates in closed form,
#
#   H(t) = lambda_0 * rho * exp(lp) * t^rho_i / rho_i,
#
# so times are still drawn exactly, with no root-finding:
#
#   T_i = ( E_i * rho_i / (lambda_0 * rho * exp(lp_i)) )^(1/rho_i),  E ~ Exp(1)
#
# At psi = 0 this reduces to the PH draw identically.  The draw requires
# rho_i > 0 for every subject, which is checked; sd(eta) is roughly 0.4-0.7 for
# the coefficient functions here, so psi up to ~0.3 is comfortable at rho ~ 1.28.
#
# Because the shape is subject-specific when psi != 0, there is NO exact AFT
# representation: the log-time error (log E)/rho_i has a subject-dependent
# scale.  Scoring such data against the psi = 0 AFT targets therefore measures
# STRUCTURAL misspecification rather than estimation error -- which is the point
# of the arm, but means beta/gamma bias should be read accordingly and the
# predictive metrics carry the comparison.
#
# RELATION TO THE AFT SCALE
# -------------------------
# Inverting the Weibull baseline gives
#   log T = -(log lambda_0)/rho - eta_Cox/rho - u_j/rho + (log E)/rho,
# so on the log-time scale beta_AFT(s) = -beta_Cox(s)/rho, the Z slopes map as
# -gamma_k/rho, the frailty variance as tau^2/rho^2, and the residual variance
# is pi^2/(6 rho^2).  `beta_scale` multiplies beta(s) before it is used: passing
# beta_scale = rho makes the implied AFT-scale curve exactly -beta(s), which is
# how a Cox arm is put on equal footing with an AFT arm specified at beta(s).
# =============================================================================

simulate_Cox <- function(
    N_total      = 4000,      # exact total sample size
    n_cluster    = 50,        # number of clusters J
    range_factor = 0.5,       # half-width of the Uniform cluster sizes, as a
                              # fraction of the mean
    tmax         = 1,         # right endpoint of the functional domain
    nS           = 401,       # grid points for s
    k0           = 20,        # Fourier basis dimension for X(s)
    alpha        = 0.7,       # eigenvalue decay: lambda_k = k^(-alpha)
    beta_type    = c("monotone", "peak1", "peak2", "wavy"),
    simulate_X   = TRUE,      # FALSE drops X(s) and the functional term entirely
    gamma        = c(0.5, 0.3, -0.2),   # c(intercept, Z1, Z2), log-hazard scale
    lambda_0     = 0.1,       # Weibull baseline rate
    rho          = 1.5,       # Weibull shape (rho = 1: exponential)
    beta_scale   = 1,         # multiplies beta(s); see "Relation to the AFT scale"
    psi          = 0,         # time-modulation of the functional effect:
                              # beta(s, t) = beta(s) * (1 + psi * log t).
                              # psi = 0 is proportional hazards; psi != 0 breaks
                              # PH.  See "Violating proportional hazards" above.
    tau          = 0.5,       # frailty SD: u_j ~ N(0, tau^2)
    u            = NULL,      # optional: supply the J cluster frailties directly
                              # instead of drawing them, e.g. to generate new
                              # subjects in KNOWN clusters.  `tau` is then unused
                              # for the draw but still recorded.
    censor_rate  = 0.25       # target censoring proportion (0 = none)
) {

  beta_type <- match.arg(beta_type)

  # ---- 0. Cluster sizes summing exactly to N_total -------------------------
  mean_nj <- N_total / n_cluster
  nj_min  <- max(1L, floor(mean_nj * (1 - range_factor)))
  nj_max  <- ceiling(mean_nj * (1 + range_factor))

  nj  <- pmax(1L, round(runif(n_cluster, nj_min, nj_max)))
  gap <- as.integer(N_total - sum(nj))
  while (gap != 0L) {
    i <- sample.int(n_cluster, 1L)
    if (gap > 0L) {
      nj[i] <- nj[i] + 1L; gap <- gap - 1L
    } else if (nj[i] > 1L) {
      nj[i] <- nj[i] - 1L; gap <- gap + 1L
    }
  }
  stopifnot(sum(nj) == N_total, all(nj >= 1L))

  N          <- sum(nj)
  cluster_id <- rep(seq_len(n_cluster), times = nj)

  # ---- 1. Scalar covariates ------------------------------------------------
  z1 <- rnorm(N, mean = 0, sd = 1)
  z2 <- rbinom(N, size = 1, prob = 0.5)

  s_grid <- seq(0, tmax, length.out = nS)

  # ---- 2. Functional predictor and its coefficient function ----------------
  #
  #   X_i(s): Karhunen-Loeve with a Fourier basis and decaying eigenvalues.
  #   beta(s): B-spline basis with the coefficient vectors from simulate_AFT.R
  #     monotone  smooth decline, no sign change
  #     peak1     single interior peak
  #     peak2     two peaks separated by a trough
  #     wavy      several sign changes -- the hardest to smooth
  # --------------------------------------------------------------------------
  wts <- trapz_weights(s_grid)

  if (simulate_X) {
    Phi_fourier      <- matrix(0, nS, k0)
    Phi_fourier[, 1] <- 1
    for (k in seq_len(floor((k0 - 1) / 2))) {
      Phi_fourier[, 2*k]     <- sqrt(2) * cos(2*pi*k*s_grid)
      Phi_fourier[, 2*k + 1] <- sqrt(2) * sin(2*pi*k*s_grid)
    }
    lambda_decay <- (seq_len(k0))^(-alpha)
    Xi           <- matrix(rnorm(N * k0), N, k0)
    sim_curves   <- Xi %*% diag(sqrt(lambda_decay)) %*% t(Phi_fourier)   # N x nS

    bs_coef <- switch(
      beta_type,
      monotone = c(0.1, -0.15, -0.35, -0.5, -0.5, -0.5),
      peak1    = c(0.0, -0.1, 0.4, 1.2, 0.5, 0.1, 0.0),
      peak2    = c(0.0, -0.1, 0.4, 0.9, 0.3, -0.15, -0.3, 0.2, 0.8, 0.3, 0.0),
      wavy     = c(0.0, 1.7, -2.1, 2.5, -0.85, 0.6, -0.4, 0.2, 0.0)
    )
    B    <- splines::bs(s_grid, df = length(bs_coef), intercept = TRUE)
    beta <- as.vector(B %*% bs_coef) * beta_scale

    num_int <- as.vector(sweep(sim_curves, 2, wts, `*`) %*% beta)
  } else {
    sim_curves <- NULL
    bs_coef    <- NULL
    beta       <- NULL
    num_int    <- numeric(N)
  }

  # ---- 3. Shared cluster frailty on the log-hazard scale -------------------
  if (is.null(u)) {
    u_j <- rnorm(n_cluster, mean = 0, sd = tau)
  } else {
    if (length(u) != n_cluster)
      stop("`u` must have length n_cluster (", n_cluster, "), got ", length(u), ".")
    u_j <- as.numeric(u)
  }
  frailty <- u_j[cluster_id]

  # ---- 4. Log-hazard linear predictor; h_0 carries the time dependence -----
  lp <- gamma[1] + gamma[2]*z1 + gamma[3]*z2 + num_int + frailty

  # ---- 5. Survival times by Weibull inverse-CDF ----------------------------
  #
  #  With beta(s, t) = beta(s)*(1 + psi*log t) the subject's shape becomes
  #  rho_i = rho + psi*eta_i and
  #    T_i = ( E_i * rho_i / (lambda_0 * rho * exp(lp_i)) )^(1/rho_i).
  #  psi = 0 gives rho_i = rho and this collapses to the PH form
  #    T_i = ( E_i / (lambda_0 * exp(lp_i)) )^(1/rho).
  rho_i <- rho + psi * num_int
  if (any(rho_i <= 0))
    stop("psi = ", psi, " makes the per-subject Weibull shape non-positive for ",
         sum(rho_i <= 0), " of ", N, " subjects (min ", round(min(rho_i), 3),
         ").  Reduce |psi| or the scale of beta(s).")

  E      <- rexp(N, rate = 1)
  T_true <- (E * rho_i / (lambda_0 * rho * exp(lp)))^(1 / rho_i)

  # ---- 6. Censoring: C ~ Uniform(0, ub) calibrated to censor_rate ----------
  if (censor_rate > 0) {
    ub <- choose_ub_unif(T_true, censor_rate)
    C  <- runif(N, 0, ub)
    Y  <- pmin(T_true, C)
  } else {
    Y <- T_true
  }
  delta <- as.integer(Y == T_true)

  # ---- 7. Assemble ---------------------------------------------------------
  sim_data <- data.frame(
    subject_id = seq_len(N),
    cluster_id = cluster_id,
    Y          = Y,
    delta      = delta,
    Z1         = z1,
    Z2         = z2,
    T_true     = T_true,
    lp         = lp,           # log hazard ratio at t = 1, excludes h_0
    rho_i      = rho_i         # per-subject Weibull shape; constant = rho if psi = 0
  )
  if (simulate_X) sim_data$X <- I(sim_curves)   # N x nS functional predictor

  list(
    data       = sim_data,
    beta       = beta,         # beta(s) on s_grid, already scaled; NULL if no X
    beta_fn    = if (simulate_X) {
                   local({ b <- beta; g <- s_grid
                           function(s) approx(g, b, xout = s, rule = 2)$y })
                 } else NULL,
    df_beta    = if (simulate_X) data.frame(s = s_grid, beta = beta) else NULL,
    s_grid     = s_grid,
    beta_type  = if (simulate_X) beta_type else NA_character_,
    beta_scale = beta_scale,
    bs_coef    = bs_coef,
    dgm        = "cox",
    gamma      = gamma,
    lambda_0   = lambda_0,
    rho        = rho,
    psi        = psi,          # 0 = proportional hazards
    tau        = tau,
    u          = u_j,
    n_subjects = nj            # realised cluster sizes
  )
}


# =============================================================================
# Helpers -- also defined in simulate_AFT.R and simulate_Cox_interaction.R;
# repeated so this file can be sourced on its own.
# =============================================================================
# Trapezoidal quadrature weights for an arbitrary grid.
trapz_weights <- function(x) {
  n  <- length(x)
  dx <- diff(x)
  w  <- numeric(n)
  w[1] <- dx[1] / 2
  w[n] <- dx[n - 1] / 2
  if (n > 2) w[2:(n - 1)] <- (dx[1:(n - 2)] + dx[2:(n - 1)]) / 2
  w
}

# Upper limit of Uniform(0, ub) censoring that hits a target censoring rate.
choose_ub_unif <- function(T_true, censor_rate, tol = 1e-8) {
  stopifnot(censor_rate > 0, censor_rate < 1)
  g  <- function(ub) mean(pmin(T_true, ub) / ub) - censor_rate
  lo <- max(tol, min(T_true[T_true > 0], na.rm = TRUE) * 1e-6)
  hi <- max(T_true, na.rm = TRUE) * 2 + 1
  if (g(lo) < 0) lo <- tol
  while (g(hi) > 0) hi <- hi * 2
  uniroot(g, lower = lo, upper = hi)$root
}


# =============================================================================
# Quick checks -- run these after changing anything above.
#
# 1. The inverse-CDF identity.  On the log scale, rho*log(T) + log(lambda_0) + lp
#    must be exactly log(E) with E ~ Exp(1), i.e. Gumbel_min(0, 1):
#    mean -gamma_EM = -0.5772, variance pi^2/6 = 1.6449.
#
# 2. Proportional hazards.  A coxph fit carrying the TRUE functional term as a
#    covariate must return gamma[2], gamma[3] and a coefficient of exactly 1 on
#    the integral.  (Fitting survreg instead recovers a Weibull shape BELOW the
#    true rho, because it has no frailty term and the unmodelled u_j inflates the
#    residual variance -- 1.495 at tau = 0, 1.03 at tau = 1.  That is a property
#    of the check, not of the simulator.)
# =============================================================================
if (FALSE) {
  library(survival)
  set.seed(99)
  s <- simulate_Cox(N_total = 40000, n_cluster = 200, nS = 101,
                    beta_type = "peak1", rho = 1.5, tau = 0.5, censor_rate = 0)

  # 1. inverse-CDF identity
  logE <- s$rho * log(s$data$T_true) + log(s$lambda_0) + s$data$lp
  cat(sprintf("mean %.4f (target %.4f) | var %.4f (target %.4f) | KS p = %.3f\n",
              mean(logE), digamma(1), var(logE), pi^2 / 6,
              ks.test(logE, function(q) 1 - exp(-exp(q)))$p.value))

  # 2. proportional hazards
  ni  <- as.vector(sweep(s$data$X, 2, trapz_weights(s$s_grid), `*`) %*% s$beta)
  fit <- coxph(Surv(T_true) ~ Z1 + Z2 + ni + frailty(cluster_id),
               data = cbind(s$data, ni = ni))
  print(data.frame(term     = c("Z1", "Z2", "int X*beta"),
                   truth    = c(s$gamma[2], s$gamma[3], 1),
                   estimate = round(as.numeric(coef(fit)[c("Z1", "Z2", "ni")]), 4)))
}
