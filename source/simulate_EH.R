####################################################################
# simulate_EH(): clustered survival data from an EXTENDED HAZARDS model
# (Etezadi-Amoli & Ciampi 1987; Chen & Jewell 2001, Biometrika 88:687-702)
# that bridges the two working models of the simulation study.
#
#   h(t | a) = h0(t * exp(kappa * a)) * exp(a),
#   a_ij     = -( z_ij' gamma + int X_ij(s) beta(s) ds + u_j ),  u_j ~ N(0, tau^2)
#
#   kappa = 1 : S(t | a) = S0(t e^{a})  ->  log T = z'gamma + int X beta + u + log T0.
#               With a log-normal baseline (log T0 ~ N(0, sigma^2)) this is EXACTLY
#               the functional frailty AFT model fitted by gibbs_functional_frailty().
#   kappa = 0 : h(t | a) = h0(t) e^{a}  ->  proportional hazards with log-hazard
#               coefficient function -beta(s), slopes -gamma[-1] and frailty -u_j,
#               i.e. EXACTLY the model fitted by the comparator fcoxFr (whose baseline
#               is unspecified).  The baseline is NOT Weibull, so there is no AFT
#               representation at kappa = 0 (a Weibull baseline would collapse the
#               whole family onto one model that is both PH and AFT).
#   0 < kappa < 1 : neither PH nor AFT; both working models are misspecified.
#
# Each working model gets its OWN nuisance component right at its own end of the
# bridge: the error distribution for the AFT model (log-normal baseline) and the
# baseline hazard for Cox (unspecified, so any h0 is "right").  A log-logistic
# baseline additionally misspecifies the AFT error at every kappa (stress test).
#
# Cumulative hazard and exact draw (no root-finding):
#   H(t | a) = exp((1 - kappa) a) * H0(t e^{kappa a})
#   T        = exp(-kappa a) * H0^{-1}( E * exp(-(1 - kappa) a) ),   E ~ Exp(1)
#
# Baselines are scaled so that SD(log T0) = sigma in both cases:
#   lognormal   : log T0 ~ N(0, sigma^2)
#   loglogistic : log T0 ~ Logistic(0, s),  s = sigma * sqrt(3) / pi
#
# Everything else (cluster sizes, Z1/Z2, Karhunen-Loeve X(s), beta(s) by
# `beta_type`, uniform censoring calibrated to `censor_rate`, the optional `u`
# argument for KNOWN-cluster test sets) mirrors simulate_AFT() / simulate_Cox(),
# so the three generators are interchangeable in the runners.
#
# Return value: list with
#   data          data.frame (Y, delta, X, W, Z1, Z2, T_true, lp, a, cluster_id, ...)
#   s_grid, beta  the grid and the TRUE beta(s) in model-(1) orientation, i.e. the
#                 exact AFT coefficient at kappa = 1; the exact log-hazard
#                 coefficient at kappa = 0 is -beta(s)
#   u             the J cluster frailties (log-time orientation)
#   surv_true     function(t, lp, u) -> matrix of TRUE conditional survival
#                 probabilities S(t | lp, u), rows = subjects, cols = t
#   surv_pop      function(t, lp, n_gh = 40) -> TRUE population-averaged survival
#                 E_u[S(t | lp, u)], u ~ N(0, tau^2), by Gauss-Hermite quadrature
#   kappa, baseline, sigma, tau, gamma, bs_coef, n_subjects, ub (censoring bound)
#
# Note: with a log-logistic baseline and kappa near 0, the UNCENSORED event-time
# distribution has a polynomial tail (T up to ~1e30 for low-risk subjects).  This
# is harmless for the analyses, which use Y = min(T, C) with C bounded, but
# survival::coxph() on uncensored T_true then fails silently: its `timefix` tie
# tolerance is relative to the time range, so most event times become ties.  Fit
# on log(T_true) or rank(T_true) if you ever need a Cox fit to uncensored times.
####################################################################

simulate_EH <- function(N_total     = 1000,
                        n_cluster   = 50,
                        range_factor = 0.5,
                        u           = NULL,
                        kappa       = 0.5,
                        baseline    = c("lognormal", "loglogistic"),
                        tmax        = 1,
                        nS          = 100,
                        k0          = 20,
                        alpha       = 0.7,
                        beta_type   = c("monotone", "peak1", "peak2", "wavy"),
                        gamma       = c(0.5, 0.3, -0.2),   # model-(1) orientation
                        sigma       = 1,                   # SD of log T0
                        tau         = 0.5,
                        censor_rate = 0.25) {

  baseline  <- match.arg(baseline)
  beta_type <- match.arg(beta_type)
  stopifnot(length(kappa) == 1L, kappa >= 0, kappa <= 1, sigma > 0, tau >= 0)

  # ---- 1. cluster sizes (identical scheme to simulate_AFT) -------------------
  mean_nj <- N_total / n_cluster
  nj_min  <- max(1L, floor(mean_nj * (1 - range_factor)))
  nj_max  <- ceiling(mean_nj * (1 + range_factor))
  nj  <- pmax(1L, round(stats::runif(n_cluster, nj_min, nj_max)))
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
  N <- N_total
  cluster_id <- rep(seq_len(n_cluster), times = nj)

  # ---- 2. scalar covariates ---------------------------------------------------
  z1 <- stats::rnorm(N)
  z2 <- stats::rbinom(N, 1, 0.5)

  # ---- 3. functional covariate (Karhunen-Loeve, Fourier; as simulate_AFT) ----
  s_grid <- seq(0, tmax, length.out = nS)
  Phi <- matrix(0, nS, k0)
  Phi[, 1] <- 1
  for (k in 1:floor((k0 - 1) / 2)) {
    Phi[, 2 * k]     <- sqrt(2) * cos(2 * pi * k * s_grid)
    Phi[, 2 * k + 1] <- sqrt(2) * sin(2 * pi * k * s_grid)
  }
  lambda_kl  <- (1:k0)^(-alpha)
  Xi         <- matrix(stats::rnorm(N * k0), N, k0)
  sim_curves <- Xi %*% diag(sqrt(lambda_kl)) %*% t(Phi)

  bs_coef <- switch(beta_type,
    monotone = c(0.1, -0.15, -0.35, -0.5, -0.5, -0.5),
    peak1    = c(0.0, -0.1, 0.4, 1.2, 0.5, 0.1, 0.0),
    peak2    = c(0.0, -0.1, 0.4, 0.9, 0.3, -0.15, -0.3, 0.2, 0.8, 0.3, 0.0),
    wavy     = c(0.0, 1.7, -2.1, 2.5, -0.85, 0.6, -0.4, 0.2, 0.0))
  B    <- splines::bs(s_grid, df = length(bs_coef), intercept = TRUE)
  beta <- as.vector(B %*% bs_coef)

  wts     <- .eh_trapz_weights(s_grid)
  Xw      <- sweep(sim_curves, 2, wts, `*`)
  num_int <- as.vector(Xw %*% beta)
  W       <- Xw %*% B

  # ---- 4. frailty (log-time orientation, as in model (1)) -------------------
  if (is.null(u)) {
    u <- stats::rnorm(n_cluster, 0, tau)
  } else {
    if (length(u) != n_cluster)
      stop("`u` must have length n_cluster (", n_cluster, "), got ", length(u), ".")
    u <- as.numeric(u)
  }

  # ---- 5. event times from the extended hazards model -----------------------
  bl <- .eh_baseline(baseline, sigma)
  lp <- gamma[1] + gamma[2] * z1 + gamma[3] * z2 + num_int   # fixed part, model (1)
  a  <- -(lp + u[cluster_id])
  E  <- stats::rexp(N)
  T_true <- exp(-kappa * a + bl$log_H0inv(log(E) - (1 - kappa) * a))
  if (any(!is.finite(T_true)) || any(T_true <= 0))
    stop("Non-finite or non-positive event times generated; check sigma/gamma/tau.")

  # ---- 6. censoring: C ~ Uniform(0, ub), ub calibrated to censor_rate -------
  if (censor_rate > 0) {
    ub <- .eh_choose_ub_unif(T_true, censor_rate)
    C  <- stats::runif(N, 0, ub)
    Y  <- pmin(T_true, C)
  } else {
    ub <- Inf
    Y  <- T_true
  }
  delta <- as.integer(T_true <= Y)

  # ---- 7. true survival functions -------------------------------------------
  # log S(t | a) = -exp( (1 - kappa) a + log H0(t e^{kappa a}) ), computed on the
  # log scale so that probabilities near 0 or 1 do not lose precision.
  surv_true <- function(t, lp, u) {
    a  <- -(lp + u)
    lt <- outer(kappa * a, log(t), `+`)               # log(t e^{kappa a})
    logH <- (1 - kappa) * a + bl$log_H0_logt(lt)       # recycled down columns
    exp(-exp(logH))
  }
  gh <- .eh_gauss_hermite(40)
  surv_pop <- function(t, lp) {
    if (tau == 0) return(surv_true(t, lp, 0))
    out <- 0
    for (k in seq_along(gh$x)) {                        # u = sqrt(2) tau x_k
      out <- out + gh$w[k] * surv_true(t, lp, sqrt(2) * tau * gh$x[k])
    }
    out / sqrt(pi)
  }

  sim_data <- data.frame(subject_id = seq_len(N),
                         cluster_id = cluster_id,
                         Y = Y, delta = delta,
                         X = I(sim_curves), W = I(W),
                         Z1 = z1, Z2 = z2,
                         T_true = T_true,
                         lp = lp,                       # fixed part, model-(1) orientation
                         a = a)                         # full EH index, incl. frailty

  list(data = sim_data, s_grid = s_grid, beta = beta, u = u,
       surv_true = surv_true, surv_pop = surv_pop,
       kappa = kappa, baseline = baseline, sigma = sigma, tau = tau,
       gamma = gamma, beta_type = beta_type, bs_coef = bs_coef,
       n_subjects = nj, ub = ub)
}

# ---- baseline distributions on the log scale ---------------------------------
# log_H0_logt(lt): log H0(t) evaluated at lt = log t
# log_H0inv(lh)  : log H0^{-1}(h) evaluated at lh = log h
.eh_baseline <- function(baseline, sigma) {
  if (baseline == "lognormal") {
    list(
      # H0(t) = -log S0(t), S0(t) = 1 - Phi(log t / sigma)
      log_H0_logt = function(lt) log(-stats::pnorm(lt / sigma, lower.tail = FALSE, log.p = TRUE)),
      # S0(t) = exp(-h)  <=>  log t = sigma * qnorm(exp(-h), lower.tail = FALSE)
      log_H0inv   = function(lh) sigma * stats::qnorm(-exp(lh), lower.tail = FALSE, log.p = TRUE)
    )
  } else {
    s <- sigma * sqrt(3) / pi
    list(
      # S0(t) = 1 / (1 + exp(log t / s))  ->  H0 = log(1 + exp(log t / s))
      log_H0_logt = function(lt) log(.eh_log1pexp(lt / s)),
      # H0^{-1}(h) = exp(s * log(exp(h) - 1))
      log_H0inv   = function(lh) s * .eh_log_expm1(exp(lh))
    )
  }
}

.eh_log1pexp  <- function(x) ifelse(x > 35, x, log1p(exp(x)))
.eh_log_expm1 <- function(h) ifelse(h > 35, h, log(expm1(h)))

.eh_trapz_weights <- function(x) {
  n <- length(x); dx <- diff(x); w <- numeric(n)
  w[1] <- dx[1] / 2; w[n] <- dx[n - 1] / 2
  if (n > 2) w[2:(n - 1)] <- (dx[1:(n - 2)] + dx[2:(n - 1)]) / 2
  w
}

.eh_choose_ub_unif <- function(T_true, censor_rate, tol = 1e-8) {
  stopifnot(censor_rate > 0, censor_rate < 1)
  g  <- function(ub) mean(pmin(T_true, ub) / ub) - censor_rate
  lo <- max(tol, min(T_true) * 1e-6)
  hi <- max(T_true) * 2 + 1
  while (g(hi) > 0) hi <- hi * 2
  stats::uniroot(g, lower = lo, upper = hi, tol = tol)$root
}

# Gauss-Hermite nodes/weights for int f(x) exp(-x^2) dx (Golub-Welsch)
.eh_gauss_hermite <- function(n) {
  i  <- seq_len(n - 1)
  J  <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i / 2)
  e  <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = sqrt(pi) * e$vectors[1, ]^2)
}
