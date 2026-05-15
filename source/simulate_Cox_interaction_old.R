# =============================================================================
# simulate_Cox_interaction.R
#
# Simulate survival data from a functional Cox model with shared Gaussian
# frailty and a bivariate coefficient surface beta(s, a), where
#   s in [0, 1]  -- the functional domain (e.g. time-of-day for MIMS)
#   a in [0, 1]  -- age scaled to the unit interval
#
# True hazard model for subject i in cluster j:
#
#   h(t_ij | u_j) = h_0(t) * exp( gamma[1] + gamma[2]*Z1_ij + gamma[3]*Z2_ij
#                                  + integral_0^1 X_ij(s) * beta(s, a_ij) ds
#                                  + u_j )
#
# where u_j ~ N(0, tau^2) is the shared log-hazard frailty.
#
# Baseline hazard h_0(t) is Weibull:
#
#   h_0(t) = lambda_0 * rho * t^(rho - 1)    [shape rho, rate lambda_0]
#   H_0(t) = lambda_0 * t^rho                 [cumulative baseline hazard]
#
# Special case rho = 1 gives an Exponential(lambda_0) baseline.
# Survival times are drawn via the inverse-CDF method:
#
#   T_ij = ( E_ij / (lambda_0 * exp(LP_ij + u_j)) )^{1/rho}
#
# where E_ij ~ Exp(1), which is equivalent to -log(Uniform(0,1)).
#
# Three coefficient surfaces are implemented, all >= 0 on [0,1]^2 and
# defined on scaled inputs (s in [0,1], a in [0,1]):
#
#   "additive"  -- beta(s,a) = s + a             (tilted plane, separable)
#   "bilinear"  -- beta(s,a) = 2*s*a             (pure interaction, separable)
#   "ridge"     -- beta(s,a) = exp(-6*(s-a)^2)   (diagonal ridge, non-separable)
#
# Helper functions trapz_weights() and choose_ub_unif() are defined at the
# bottom of this file (identical to those in simulate_AFT_interaction.R).
# =============================================================================
simulate_Cox_interaction <- function(
    n_cluster   = 100,
    n_subject   = 5,
    tmax        = 1,          # right endpoint of functional domain s in [0, tmax]
    nS          = 401,        # number of grid points for s
    k0          = 20,         # Fourier basis dimension for simulating X(s)
    alpha       = 0.7,        # eigenvalue decay rate: lambda_k = k^{-alpha}
    beta_type   = c("additive", "bilinear", "ridge"),
    age_range   = c(40, 85),  # raw age drawn from Uniform(age_range[1], age_range[2])
    gamma       = c(0.5, 0.3, -0.2), # c(intercept, Z1 coef, Z2 coef)
    lambda_0    = 0.1,        # Weibull baseline rate  (h_0(t) = lambda_0*rho*t^{rho-1})
    rho         = 1.5,        # Weibull shape  (rho=1 -> exponential baseline)
    tau         = 0.5,        # frailty SD  (u_j ~ N(0, tau^2))
    censor_rate = 0.25        # target censoring proportion (0 = no censoring)
) {

  beta_type <- match.arg(beta_type)

  # ---------------------------------------------------------------------------
  # 0.  Dimensions and indexing
  # ---------------------------------------------------------------------------
  N          <- n_cluster * n_subject
  cluster_id <- rep(seq_len(n_cluster), each = n_subject)

  # ---------------------------------------------------------------------------
  # 1.  Scalar covariates
  # ---------------------------------------------------------------------------
  z1 <- rnorm(N, mean = 0, sd = 1)
  z2 <- rbinom(N, size = 1, prob = 0.5)

  # ---------------------------------------------------------------------------
  # 2.  Age: simulate raw, then scale to a in [0, 1] for surface evaluation
  # ---------------------------------------------------------------------------
  age_raw    <- runif(N, min = age_range[1], max = age_range[2])
  age_scaled <- (age_raw - age_range[1]) / diff(age_range)
  # age_scaled = 0  <->  youngest (age_range[1])
  # age_scaled = 1  <->  oldest   (age_range[2])

  # ---------------------------------------------------------------------------
  # 3.  Functional predictor X_i(s) via Fourier basis + decaying eigenvalues
  # ---------------------------------------------------------------------------
  s_grid <- seq(0, tmax, length.out = nS)

  # Fourier basis matrix (nS x k0)
  Phi_fourier      <- matrix(0, nS, k0)
  Phi_fourier[, 1] <- 1
  for (k in seq_len(floor((k0 - 1) / 2))) {
    Phi_fourier[, 2*k]     <- sqrt(2) * cos(2*pi*k*s_grid)
    Phi_fourier[, 2*k + 1] <- sqrt(2) * sin(2*pi*k*s_grid)
  }

  lambda_decay <- (seq_len(k0))^(-alpha)          # eigenvalue decay
  Xi           <- matrix(rnorm(N * k0), N, k0)    # random scores
  sim_curves   <- Xi %*% diag(sqrt(lambda_decay)) %*% t(Phi_fourier)  # N x nS

  # ---------------------------------------------------------------------------
  # 4.  True bivariate coefficient surface beta(s, a)
  #
  #   s in [0, 1] -- normalised functional domain
  #   a in [0, 1] -- normalised age
  #   All surfaces are >= 0 everywhere and contain no additive constant,
  #   so absorb.cons in the fitting function does not remove signal.
  # ---------------------------------------------------------------------------

  if (beta_type == "additive") {

    # ------------------------------------------------------------------
    # Surface 1: additive plane  beta(s, a) = s + a
    #
    #   Separable, no interaction. The s and a effects are independent.
    #   Lies exactly in the span of the marginal basis functions --
    #   the easiest case for the tensor smooth to represent.
    #
    #   Symmetric: beta(s,a) = beta(a,s)
    #   Range: [0, 2].  Zero only at the corner (s=0, a=0).
    #   Domain mean: E[s] + E[a] = 1.
    # ------------------------------------------------------------------
    beta_fn <- function(s, a) s + a

  } else if (beta_type == "bilinear") {

    # ------------------------------------------------------------------
    # Surface 2: bilinear interaction  beta(s, a) = 2*s*a
    #
    #   Non-zero only when BOTH s and a are large. Requires cross-product
    #   tensor basis functions; cannot be represented by marginals alone.
    #
    #   Symmetric: beta(s,a) = beta(a,s)
    #   Range: [0, 2].  Zero along s=0 or a=0.
    #   Domain mean: 2*E[s]*E[a] = 0.5.
    # ------------------------------------------------------------------
    beta_fn <- function(s, a) 2 * s * a

  } else if (beta_type == "ridge") {

    # ------------------------------------------------------------------
    # Surface 3: Gaussian diagonal ridge  beta(s, a) = exp(-6*(s-a)^2)
    #
    #   Large along the main diagonal s = a; decays rapidly off-diagonal.
    #   Genuinely non-separable: no f(s)*g(a) factorisation exists.
    #
    #   Symmetric: beta(s,a) = beta(a,s)  [exact]
    #   Range: [exp(-6), 1] ~= [0.002, 1].  Strictly positive everywhere.
    #   Domain mean: ~0.56.
    # ------------------------------------------------------------------
    beta_fn <- function(s, a) exp(-6 * (s - a)^2)

  }

  # ---------------------------------------------------------------------------
  # 5.  Numerical integration:  num_int_i = sum_s X_i(s) * beta(s, a_i) * w(s)
  #
  #   Because beta depends on each subject's age, the surface is evaluated
  #   separately for each i.
  # ---------------------------------------------------------------------------
  wts     <- trapz_weights(s_grid)
  num_int <- numeric(N)

  for (i in seq_len(N)) {
    beta_i     <- beta_fn(s_grid, age_scaled[i])
    num_int[i] <- sum(sim_curves[i, ] * beta_i * wts)
  }

  # ---------------------------------------------------------------------------
  # 6.  Shared cluster frailty  u_j ~ N(0, tau^2)  on the log-hazard scale
  # ---------------------------------------------------------------------------
  u_j     <- rnorm(n_cluster, mean = 0, sd = tau)
  frailty <- u_j[cluster_id]

  # ---------------------------------------------------------------------------
  # 7.  Log-hazard linear predictor (same for all t; h_0 carries the time dep.)
  # ---------------------------------------------------------------------------
  lp <- gamma[1] + gamma[2]*z1 + gamma[3]*z2 + num_int + frailty
  # lp = log{ h(t) / h_0(t) }  for each subject

  # ---------------------------------------------------------------------------
  # 8.  Survival times via Weibull inverse-CDF
  #
  #   h(t) = lambda_0 * rho * t^{rho-1} * exp(lp)
  #   H(t) = lambda_0 * t^rho * exp(lp)
  #   S(t) = exp(-H(t))
  #
  #   Inverse CDF:  T = [ -log(U) / (lambda_0 * exp(lp)) ]^{1/rho}
  #                   = [ E / (lambda_0 * exp(lp)) ]^{1/rho}
  #   where E ~ Exp(1) = -log(Uniform(0,1)).
  # ---------------------------------------------------------------------------
  E      <- rexp(N, rate = 1)
  T_true <- (E / (lambda_0 * exp(lp)))^(1 / rho)

  # ---------------------------------------------------------------------------
  # 9.  Censoring via Uniform(0, ub) calibrated to hit censor_rate
  # ---------------------------------------------------------------------------
  if (censor_rate > 0) {
    ub <- choose_ub_unif(T_true, censor_rate)
    C  <- runif(N, 0, ub)
    Y  <- pmin(T_true, C)
  } else {
    Y <- T_true
  }
  delta <- as.integer(Y == T_true)

  # ---------------------------------------------------------------------------
  # 10.  Data frame for model fitting
  #
  #   The integrated tensor basis W is NOT precomputed here because it is
  #   subject-specific. Pass X and age to the fitting function directly.
  # ---------------------------------------------------------------------------
  sim_data <- data.frame(
    subject_id = seq_len(N),
    cluster_id = cluster_id,
    Y          = Y,
    delta      = delta,
    X          = I(sim_curves),   # N x nS functional predictor matrix
    Z1         = z1,
    Z2         = z2,
    age        = age_raw,         # raw age (for the fitting function)
    age_scaled = age_scaled,      # a in [0,1]  (for DGM evaluation)
    T_true     = T_true,
    lp         = lp               # log-hazard ratio (excludes h_0)
  )

  # ---------------------------------------------------------------------------
  # 11.  True surface on a plotting grid  (nS x 51)
  # ---------------------------------------------------------------------------
  a_eval   <- seq(0, 1, length.out = 51)
  age_eval <- a_eval * diff(age_range) + age_range[1]

  beta_surface_mat <- outer(s_grid, a_eval, FUN = beta_fn)   # nS x 51

  beta_surface <- data.frame(
    s          = rep(s_grid,   times = length(a_eval)),
    age_scaled = rep(a_eval,   each  = nS),
    age        = rep(age_eval, each  = nS),
    beta       = as.vector(beta_surface_mat)
  )

  # ---------------------------------------------------------------------------
  # 12.  Return
  # ---------------------------------------------------------------------------
  list(
    data             = sim_data,
    beta_fn          = beta_fn,          # closure: beta_fn(s, a) for any s, a
    beta_surface     = beta_surface,     # long data frame for ggplot
    beta_surface_mat = beta_surface_mat, # nS x 51 matrix
    s_grid           = s_grid,
    beta_type        = beta_type,
    dgm              = "cox",
    age_range        = age_range,
    gamma            = gamma,
    lambda_0         = lambda_0,
    rho              = rho,
    tau              = tau,
    u                = u_j
  )
}


# =============================================================================
# Helper: trapezoidal quadrature weights for an arbitrary grid x
# =============================================================================
trapz_weights <- function(x) {
  n  <- length(x)
  dx <- diff(x)
  w  <- numeric(n)
  w[1] <- dx[1] / 2
  w[n] <- dx[n - 1] / 2
  if (n > 2) w[2:(n - 1)] <- (dx[1:(n - 2)] + dx[2:(n - 1)]) / 2
  w
}


# =============================================================================
# Helper: find the upper limit ub of Uniform(0, ub) censoring so that
#         E[1(T <= C)] = censor_rate  (solved numerically via uniroot)
# =============================================================================
choose_ub_unif <- function(T_true, censor_rate, tol = 1e-8) {
  stopifnot(censor_rate > 0, censor_rate < 1)

  g <- function(ub) mean(pmin(T_true, ub) / ub) - censor_rate

  lo <- max(tol, min(T_true[T_true > 0], na.rm = TRUE) * 1e-6)
  hi <- max(T_true, na.rm = TRUE) * 2 + 1

  if (g(lo) < 0) lo <- tol
  while (g(hi) > 0) hi <- hi * 2

  uniroot(g, lower = lo, upper = hi)$root
}
