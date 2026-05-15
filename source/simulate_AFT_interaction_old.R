# =============================================================================
# simulate_AFT_interaction.R
#
# Simulate survival data from a lognormal (or log-logistic) AFT model whose
# true functional effect is a BIVARIATE coefficient surface beta(s, a), where
#   s in [0, 1]  -- the functional domain (e.g. time-of-day for MIMS)
#   a in [0, 1]  -- age scaled to the unit interval
#
# True linear predictor for subject i in cluster j:
#
#   log T_ij = gamma[1] + gamma[2]*z1_ij + gamma[3]*z2_ij
#              + integral_0^1 X_ij(s) * beta(s, a_ij) ds
#              + u_j  +  sigma * epsilon_ij
#
# where u_j ~ N(0, tau^2) is the shared cluster frailty.
#
# Two surface types are implemented, both >= 0 everywhere and defined on the
# scaled inputs (s in [0,1], a in [0,1]):
#
#   "additive"  -- beta(s,a) = 1 + s + a   (tilted plane, no interaction term)
#   "bilinear"  -- beta(s,a) = 1 + 2*s*a   (pure s x a interaction)
#
# Helper functions trapz_weights() and choose_ub_unif() are appended at the
# bottom of this file.
# =============================================================================

simulate_AFT_interaction <- function(
    family       = c("lognormal", "loglogistic"),
    n_cluster    = 100,
    n_subject    = 5,
    tmax         = 1,         # right endpoint of functional domain s in [0, tmax]
    nS           = 401,       # number of grid points for s
    k0           = 20,        # Fourier basis dimension for simulating X(s)
    alpha        = 0.7,       # eigenvalue decay rate: lambda_k = k^{-alpha}
    beta_type = c("additive", "bilinear", "ridge"),
    age_range    = c(40, 85), # raw age drawn from Uniform(age_range[1], age_range[2])
    gamma        = c(0.5, 0.3, -0.2), # c(intercept, z1 coef, z2 coef)
    sigma        = 0.2,       # error SD
    tau          = 0.5,       # frailty SD
    censor_rate  = 0.25       # target censoring proportion (0 = no censoring)
) {
  
  family       <- match.arg(family)
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
  
  lambda_decay <- (seq_len(k0))^(-alpha)         # eigenvalue decay
  Xi           <- matrix(rnorm(N * k0), N, k0)   # random scores
  sim_curves   <- Xi %*% diag(sqrt(lambda_decay)) %*% t(Phi_fourier)  # N x nS
  
  # ---------------------------------------------------------------------------
  # 4.  True bivariate coefficient surface  beta(s, a)
  #
  #  Both surfaces use:
  #    s in [0, 1]  (normalised functional domain)
  #    a in [0, 1]  (normalised age)
  #  and are >= 0 everywhere by construction.
  # ---------------------------------------------------------------------------
  
  if (beta_type == "additive") {
    
    # ------------------------------------------------------------------
    # Surface 1: additive plane  beta(s, a) = s + a
    #
    #   Separable, no interaction. The effect of s does not depend on a
    #   and vice versa. Lies exactly in the span of the marginal basis
    #   functions -- the easiest case for the tensor model.
    #
    #   beta(s, a) = beta(a, s)  [symmetric in s and a]
    #   Range : [0, 2].  Zero only at the corner (s=0, a=0).
    #   Domain mean : E[s] + E[a] = 0.5 + 0.5 = 1.
    # ------------------------------------------------------------------
    beta_fn <- function(s, a) s + a
    
  } else if (beta_type == "bilinear") {
    
    # ------------------------------------------------------------------
    # Surface 2: bilinear interaction  beta(s, a) = 2*s*a
    #
    #   Non-separable multiplicative interaction: the surface is large
    #   only when BOTH s and a are large, and zero along both axes.
    #   Requires the cross-product tensor basis functions to represent.
    #
    #   beta(s, a) = beta(a, s)  [symmetric in s and a]
    #   Range : [0, 2].  Zero along s=0 or a=0.
    #   Domain mean : 2 * E[s] * E[a] = 0.5.
    # ------------------------------------------------------------------
    beta_fn <- function(s, a) 2 * s * a
    
  } else if (beta_type == "ridge") {
    
    # ------------------------------------------------------------------
    # Surface 3: Gaussian diagonal ridge  beta(s, a) = exp(-6*(s-a)^2)
    #
    #   The surface is concentrated along the main diagonal s = a, where
    #   beta = 1, and decays rapidly away from it (half-max at |s-a|=0.34).
    #   Unlike additive and bilinear, this surface is high at BOTH corners
    #   (0,0) and (1,1) and low at the off-diagonal corners (1,0) and (0,1).
    #   It is genuinely non-separable: no f(s)*g(a) factorisation exists.
    #
    #   beta(s, a) = beta(a, s)  [symmetric in s and a, exact]
    #   Range : [exp(-6), 1] ≈ [0.002, 1].  Strictly positive everywhere.
    #   Domain mean : ~0.56  (from Monte Carlo).
    # ------------------------------------------------------------------
    beta_fn <- function(s, a) exp(-6 * (s - a)^2)
    
  }
  
  # ---------------------------------------------------------------------------
  # 5.  Numerical integration:  int_i = sum_s X_i(s) * beta(s, a_i) * w(s)
  #
  #  Because beta depends on the individual's age, the surface must be
  #  evaluated separately for each subject i.
  # ---------------------------------------------------------------------------
  wts     <- trapz_weights(s_grid)   # trapezoidal quadrature weights, length nS
  num_int <- numeric(N)
  
  for (i in seq_len(N)) {
    beta_i   <- beta_fn(s_grid, age_scaled[i])    # evaluate surface at age_i
    num_int[i] <- sum(sim_curves[i, ] * beta_i * wts)
  }
  
  # ---------------------------------------------------------------------------
  # 6.  Error term
  # ---------------------------------------------------------------------------
  z_err <- if (family == "lognormal") rnorm(N) else stats::rlogis(N)
  
  # ---------------------------------------------------------------------------
  # 7.  Shared cluster frailty
  # ---------------------------------------------------------------------------
  u_j    <- rnorm(n_cluster, mean = 0, sd = tau)
  frailty <- u_j[cluster_id]
  
  # ---------------------------------------------------------------------------
  # 8.  Linear predictor and true survival times
  # ---------------------------------------------------------------------------
  lp     <- gamma[1] + gamma[2]*z1 + gamma[3]*z2 + num_int
  T_true <- exp(lp + frailty + sigma * z_err)
  
  # ---------------------------------------------------------------------------
  # 9.  Censoring via Uniform(0, ub) where ub is chosen to hit censor_rate
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
  #  Note: W (the integrated tensor basis matrix) is NOT precomputed here
  #  because it is subject-specific in the interaction model; gibbs_frailty_
  #  interaction() constructs it internally using PredictMat().
  #  Pass sim_data$X and sim_data$age directly to that function.
  # ---------------------------------------------------------------------------
  sim_data <- data.frame(
    subject_id = seq_len(N),
    cluster_id = cluster_id,
    Y          = Y,
    delta      = delta,
    X          = I(sim_curves),   # N x nS functional predictor matrix
    Z1         = z1,
    Z2         = z2,
    age        = age_raw,         # raw age (pass to gibbs_frailty_interaction)
    age_scaled = age_scaled,      # a in [0,1]  (used for DGM / diagnostics)
    T_true     = T_true,
    lp         = lp
  )
  
  # ---------------------------------------------------------------------------
  # 11.  True surface on a plotting grid  (nS x n_age_eval)
  #
  #  beta_surface$beta  -- vectorised values of beta(s, a)
  #  beta_surface_mat   -- matrix form (nS rows = s, n_age_eval cols = age)
  # ---------------------------------------------------------------------------
  a_eval   <- seq(0, 1, length.out = 51)
  age_eval <- a_eval * diff(age_range) + age_range[1]
  
  beta_surface_mat <- outer(s_grid, a_eval, FUN = beta_fn)  # nS x 51
  
  beta_surface <- data.frame(
    s          = rep(s_grid,  times = length(a_eval)),
    age_scaled = rep(a_eval,  each  = nS),
    age        = rep(age_eval, each = nS),
    beta       = as.vector(beta_surface_mat)
  )
  
  # ---------------------------------------------------------------------------
  # 12.  Return
  # ---------------------------------------------------------------------------
  list(
    data             = sim_data,
    beta_fn          = beta_fn,         # closure: beta_fn(s, a) for any s, a
    beta_surface     = beta_surface,    # long data frame for ggplot
    beta_surface_mat = beta_surface_mat,# nS x 51 matrix form
    s_grid           = s_grid,
    beta_type     = beta_type,
    family           = family,
    age_range        = age_range,
    gamma            = gamma,
    sigma            = sigma,
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