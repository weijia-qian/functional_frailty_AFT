####################################################################
# Weijia Qian
# Jan 10, 2026
#
# Simulate survival data under a functional frailty AFT model with a
# NON-INTERACTION coefficient function beta(s), fit gibbs_functional_frailty(),
# and extract estimation + predictive metrics.
#
# This is the beta(s) counterpart of main_simulations_interaction.R and uses the
# same design grid, test sets, targets and metrics, so the two studies are
# directly comparable.  The differences are exactly two: the DGM is
# simulate_AFT() (beta depends on s only, and no age is simulated), and the
# fitted model is gibbs_functional_frailty() rather than the tensor version.
#
# The "cox_tv" arm is the same Weibull PH model with the functional effect
# drifting in time, beta(s, t) = beta(s)*(1 + psi*log t), which VIOLATES
# proportional hazards: each subject's Weibull shape becomes rho + psi*eta.  It
# is the only arm whose truth has no AFT representation at all, so its beta /
# gamma / tau2 / sigma2 targets are the psi = 0 ones and measure STRUCTURAL
# misspecification rather than estimation error.  Read the predictive metrics
# (C-index, IBS) as the meaningful comparison for that arm.
#
# Everything else is scored on the AFT (log-time) scale.  Both families here are AFT
# DGMs, so the DGM value and the target coincide -- unlike the interaction study,
# there is no Cox arm and hence no log-hazard -> log-time transform.
#
# Per-iteration output (result$coef)
# ----------------------------------
#   beta_{ise,ise_rel,bias_mean,bias_max}   ise_rel = ise / mean(beta_true^2)
#   beta_cover_pw / beta_cover_sim          pointwise / simultaneous 95% coverage
#   beta_cma_qd                             simultaneous critical value
#   gamma_{true,est,bias,se,cover}_k
#   tau2_{true,est,bias,se,cover}
#   sigma2_{true,est,bias,se,cover}         target is the log-time residual variance
#   c_index, ibs                            marginal: new clusters, frailty = 0
#   {c_index,ibs}_known                     conditional: known clusters, frailty = u_hat
#   {c_index,ibs}_known_marg                same subjects, frailty = 0
#   cor_u                                   cor(u_hat_j, u_j)
#   risk_sd, risk_n_unique                  guards against a collapsed risk score
#   lf_*                                    linear functional Cox + frailty
#                                           comparator (fcoxFr): predictive
#                                           metrics, the frailty variance, gamma,
#                                           and the SHAPE of beta_hat.  It has no
#                                           interval estimates, so no coverage.
#
# result$info carries the design, event rate and sampler seconds; result$beta
# the fitted beta(s) on the evaluation grid.
####################################################################

suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survAUC))    # UnoC (IPCW C-index)
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(coxme))      # comparator: frailty Cox
suppressPackageStartupMessages(library(fda))        # comparator: FPCA basis

wd = getwd()

if(substring(wd, 2, 6) == "Users"){
  doLocal = TRUE
}else{
  doLocal = FALSE
}

###############################################################
## define or source functions used in code below
###############################################################
source(here("source", "gibbs.R"))               # gibbs_functional_frailty, cma_band
source(here("source", "predict_gibbs.R"))
source(here("source", "cal_Brier.R"))
source(here("source", "simulate_AFT.R"))
source(here("source", "simulate_Cox.R"))
source(here("source", "fcoxFr_supp.R"))   # comparator: fpcr / fpcr_predict

###############################################################
## set simulation design elements
###############################################################
family = c("lognormal", "cox", "cox_tv")
# family = c("lognormal", "loglogistic", "cox", "cox_tv")
N_total = c(4000)
n_cluster = c(25, 50, 75)
nS = c(100)
beta_type = c("peak1")
# beta_type = c("monotone", "peak1", "peak2", "wavy")
tau = c(0.2, 0.5, 1)  # AFT-scale (log-time) frailty SD target, common to both families
sigma = c(1)       # AFT-scale (log-time) RESIDUAL SD target, common to both
                   # families.  Each DGM takes a different scale parameter to
                   # hit it.  See below.
censor_rate = c(0.5)

# Strength of the PH violation in the "cox_tv" arm: beta(s, t) = beta(s)*(1 + psi*log t),
# so the subject's Weibull shape becomes rho + psi*eta.  psi = 0 would BE the
# "cox" arm; 0.3 gives a violation cox.zph rejects at p < 1e-40 (N = 8000, no
# frailty) while keeping every rho_i comfortably positive.
PSI_TV = 0.3
N_iter = 500
N_test = 1000      # out-of-sample test set size (new subjects, NEW clusters)
N_known = 1000     # out-of-sample test set size (new subjects, KNOWN clusters);
                   # drawn per iteration sharing the training data's frailties

# Fitting configuration.  K_fit / bs_fit match the s-margin of the interaction
# study's tensor, so the two studies use the same basis in s.
K_fit      <- 15
bs_fit     <- "cr"
n_iter_fit <- 20000
n_burn_fit <- 10000

params = expand.grid(family = family,
                     N_total = N_total,
                     n_cluster = n_cluster,
                     nS = nS,
                     beta_type = beta_type,
                     tau = tau,
                     sigma = sigma,
                     censor_rate = censor_rate)

## define number of simulations and parameter scenarios
if(doLocal) {
  scenario = 1
  N_iter = 1
}else{
  # defined from batch script params
  scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
}

# Validate here rather than letting an out-of-range or NA index propagate:
# params$family[scenario] would silently be NA and the first `if` that touches
# `family` fails with "missing value where TRUE/FALSE needed", far from the
# actual mistake.  nrow(params) changes whenever the design grid changes.
if (length(scenario) != 1L || is.na(scenario) ||
    scenario < 1L || scenario > nrow(params)) {
  stop(sprintf(
    "scenario must be a single integer in 1..%d, got: %s\n  (doLocal = %s, commandArgs = %s)",
    nrow(params), paste(deparse(scenario), collapse = ""), doLocal,
    paste(deparse(commandArgs(trailingOnly = TRUE)), collapse = "")))
}

family      = params$family[scenario]
N_total     = params$N_total[scenario]
n_cluster   = params$n_cluster[scenario]
nS          = params$nS[scenario]
beta_type   = params$beta_type[scenario]
tau         = params$tau[scenario]
sigma       = params$sigma[scenario]
censor_rate = params$censor_rate[scenario]

# Both Cox arms share the log-hazard scale and therefore all the scale
# corrections; they differ only in whether the effect drifts with time.
is_cox <- as.character(family) %in% c("cox", "cox_tv")

# sigma is the shared log-time RESIDUAL SD target, so both families are compared
# at matched signal-to-noise and only the error SHAPE differs.  Each DGM needs a
# different scale parameter to hit it, since sigma_err = sigma_dgm * sd(epsilon):
#   lognormal    sd(N(0,1))        = 1          -> sigma_dgm = sigma
#   loglogistic  sd(Logistic(0,1)) = pi/sqrt(3) -> sigma_dgm = sigma*sqrt(3)/pi
#   cox          the log-time error is Gumbel_min(0,1)/rho, so rho is not free
sigma_dgm <- if (is_cox) NA_real_ else switch(as.character(family),
                    lognormal   = sigma,
                    loglogistic = sigma * sqrt(3) / pi)

rho_dgm <- if (is_cox) pi / (sqrt(6) * sigma) else NA_real_

# Guard: the derivation above must leave every family with the SAME log-time
# residual variance, sigma^2.  If this fires, sigma_dgm / rho_dgm and the
# sigma2_true switch inside the loop have drifted out of sync.
# For cox_tv the shape is subject-specific (rho + psi*eta), so pi^2/(6 rho^2) is
# the psi = 0 value, i.e. a reference point rather than an exact residual
# variance -- see the note on structural misspecification below.
stopifnot(isTRUE(all.equal(
  if (is_cox) pi^2 / (6 * rho_dgm^2) else switch(as.character(family),
         lognormal   = sigma_dgm^2,
         loglogistic = sigma_dgm^2 * pi^2 / 3),
  sigma^2)))

# gamma is specified on the AFT scale.  gamma_AFT = -gamma_Cox/rho, so the Cox
# arm is handed rho*gamma to give its Z slopes the same AFT-scale magnitude.
# Only the slopes are matched; the intercept is a pure time-scale offset that
# the censoring calibration and the Brier t-grid both absorb.
gamma_aft <- c(0.5, 0.3, -0.2)      # intercept, Z1, Z2 -- AFT (log-time) scale
gamma_dgm <- if (is_cox) rho_dgm * gamma_aft else gamma_aft

# tau is the shared AFT-scale frailty SD target.  u_AFT = -u_Cox/rho, so the Cox
# arm is simulated at rho_dgm * tau to give an implied AFT-scale SD of exactly tau.
tau_dgm <- if (is_cox) rho_dgm * tau else tau

# AFT-scale targets under the Cox DGM.  Inverting the Weibull baseline gives
#   log T = -(log lambda_0)/rho - eta_Cox/rho - u/rho + (log E)/rho,  E ~ Exp(1)
# so every effect maps to log time by division by -rho.  The intercept also
# absorbs the baseline scale -(log lambda_0)/rho and, because the fitted
# Gaussian error is mean zero while (log E)/rho has mean -gamma_EM/rho, that
# offset too.
EULER_MASCHERONI <- -digamma(1)          # 0.5772157
LAMBDA_0         <- 0.1                  # Weibull baseline rate for the Cox DGM

aft_beta_target <- function(beta_cox, rho) -beta_cox / rho

aft_gamma_target <- function(gamma_cox, rho, lambda_0) {
  c(-(log(lambda_0) + gamma_cox[1] + EULER_MASCHERONI) / rho,
    -gamma_cox[-1] / rho)
}

# One DGM call for whichever family this scenario is.  `u` and N are the only
# things that vary between the training, marginal-test and known-cluster draws.
sim_dgm <- function(N, J, u = NULL) {
  if (is_cox) {
    simulate_Cox(N_total = N, n_cluster = J, nS = nS,
                 beta_type = as.character(beta_type), gamma = gamma_dgm,
                 lambda_0 = LAMBDA_0, rho = rho_dgm,
                 beta_scale = rho_dgm,   # implied AFT curve = -beta(s)
                 psi = if (as.character(family) == "cox_tv") PSI_TV else 0,
                 tau = tau_dgm, u = u, censor_rate = censor_rate)
  } else {
    simulate_AFT(family = as.character(family), N_total = N, n_cluster = J,
                 nS = nS, beta_type = as.character(beta_type), gamma = gamma_dgm,
                 sigma = sigma_dgm, tau = tau_dgm, u = u, censor_rate = censor_rate)
  }
}

# simulate_AFT() returns beta on s_grid via $coefficients$beta; simulate_Cox()
# via $beta.  This hides the difference.
dgm_beta   <- function(s) if (is_cox) s$beta else s$coefficients$beta
dgm_s_grid <- function(s) s$s_grid

###############################################################
## run simulations
###############################################################

coef_list <- vector("list", length = N_iter)
info_list <- vector("list", length = N_iter)

###############################################################
## Marginal test set: one draw shared by every iteration.
## New clusters -> frailty = 0 at prediction time.
##
## The reseed is NOT for reproducibility -- set.seed(scenario) already gives
## that.  It decouples the test set from N_iter: without it the RNG state here
## would depend on how many draws sample.int(1e8, N_iter) consumed.
###############################################################
set.seed(scenario)

# Deliberately NOT a mirror of the training design.  Predictions here set the
# frailty to 0, so the test clusters exist only to inject between-cluster
# variance; since this set is drawn ONCE and reused by every iteration, few
# clusters would mean a small sample of u_j whose idiosyncrasy biases the whole
# scenario and never averages out.
n_cluster_test <- max(50L, N_test %/% 5L)

sim_test  <- sim_dgm(N_test, n_cluster_test)
test_data <- sim_test$data
Z_test    <- model.matrix(~ Z1 + Z2, data = test_data)
s_grid    <- dgm_s_grid(sim_test)
cat("Test set generated: N =", nrow(test_data),
    "| event rate =", round(mean(test_data$delta), 3), "\n")

# tgrid: 5th-95th percentiles of test EVENT times (not censored), so the grid
# spans the event-time distribution including the tails.  Fixed across iterations.
event_times_test <- test_data$Y[test_data$delta == 1]
tgrid <- quantile(event_times_test, probs = seq(0.05, 0.95, by = 0.05))
tgrid <- sort(unique(tgrid))

# Trim tgrid to where the IPCW weight 1/G(t) is estimable, using the test set's
# own censoring distribution.  Both are fixed across iterations, so the horizon
# is constant within a scenario and every replicate estimates the same
# functional.  Estimating G on the training data instead is what produced the
# IBS blow-ups: its reverse KM hits 0 at max(time_train).
eps_G      <- 0.05
km_cens_te <- survival::survfit(survival::Surv(test_data$Y, 1 - test_data$delta) ~ 1)
G_test     <- stats::stepfun(km_cens_te$time, c(1, km_cens_te$surv))
tgrid      <- tgrid[G_test(tgrid) >= eps_G]
stopifnot(length(tgrid) >= 2L)

# Shared truncation time for both predictive metrics.  Uno's C-index weights by
# 1/G(T)^2, which is least stable at the very end of follow-up, so pass an
# explicit tau rather than letting UnoC default to max(test time).
tau_eval <- max(tgrid)
cat("Evaluation horizon: tgrid has", length(tgrid), "points, tau =",
    round(tau_eval, 4), "| G(tau) =", round(G_test(tau_eval), 4), "\n")

Surv_test <- survival::Surv(test_data$Y, test_data$delta)

# beta(s) is evaluated on the simulation grid itself, which is constant within a
# scenario, as is the truth.
# AFT-scale target curve: unchanged for the AFT families, -beta_Cox/rho for the
# Cox DGM, whose beta lives on the log-hazard scale.  Because the Cox arm is
# specified at rho*beta, this evaluates to -beta(s) -- the same magnitude as the
# AFT arms, sign-flipped as the transform requires.
beta_true_dgm <- dgm_beta(sim_test)
beta_true <- if (is_cox) {
  aft_beta_target(beta_true_dgm, rho_dgm)
} else {
  beta_true_dgm
}
n_point   <- length(s_grid)

mat_grid <- function() matrix(NA_real_, N_iter, n_point)
bmat     <- list(mean = mat_grid(), sd = mat_grid(),
                 q025 = mat_grid(), q975 = mat_grid())

# Per-iteration seeds
seeds <- sample.int(1e8, N_iter)

for(iter in 1:N_iter){
  cat("Scenario", scenario, "| Iteration", iter, "\n")
  set.seed(seeds[[iter]])

  # simulate training data -- retry up to 10 times if cluster-size sampling fails
  sim_data <- NULL
  for (.attempt in seq_len(10L)) {
    sim_data <- tryCatch(
      sim_dgm(N_total, n_cluster),
      error = function(e) {
        warning(sprintf("Simulation attempt %d failed (iter %d): %s",
                        .attempt, iter, e$message))
        NULL
      })
    if (!is.null(sim_data)) break
  }
  if (is.null(sim_data)) {
    warning(sprintf("Iter %d: all simulation attempts failed; skipping.", iter))
    next
  }

  res <- tryCatch({

    ###############################################################
    ## fit the functional frailty AFT model
    ###############################################################
    data <- sim_data$data
    Z    <- model.matrix(~ Z1 + Z2, data = data)

    tic()
    fit <- gibbs_functional_frailty(
      time = data$Y, status = data$delta, cluster_id = data$cluster_id,
      Z = Z, X = data$X, s_grid = s_grid,
      # tuning / priors -- matched to the interaction study
      K = K_fit, basis_type = bs_fit,
      lambda_init = 1, A_lambda = 1, B_lambda = 1,
      var_gamma = 100,
      # IG(0.01, 0.01) for tau^2, not IG(3, 2).  IG(3, 2) has prior mean 1, which
      # swamps the frailty likelihood when tau is small: at tau = 0.2 it gave
      # tau2_hat ~ 0.15-0.21 against a truth of 0.04 and coverage of exactly 0,
      # at both J = 25 and J = 50.  The near-flat prior recovers 0.043-0.060 with
      # coverage 0.92-1.00 and costs nothing at tau = 0.5 or 1 (coverage 1.00
      # either way at 30 reps).  coxme recovers 0.039 on the same data, so the
      # information is in the likelihood and only the prior was obstructing it.
      # Note IG(3, 2)'s prior mean of 1 also happens to equal the tau = 1 truth,
      # which flattered that cell.
      A_tau2 = 0.01, B_tau2 = 0.01,
      A_sigma2 = 3, B_sigma2 = 2,  # IG(A,B) for sigma^2
      n_iter = n_iter_fit, n_burn = n_burn_fit, n_thin = 1, verbose = TRUE
    )
    time_stamp <- toc(quiet = TRUE)
    time <- time_stamp$toc - time_stamp$tic

    ###############################################################
    ## beta(s) estimation metrics
    ###############################################################
    # gibbs_functional_frailty() already returns the posterior summaries and the
    # simultaneous band on s_grid.
    beta_mean <- fit$beta_mean
    beta_sd   <- fit$beta_cma_sd
    beta_q025 <- fit$beta_q025
    beta_q975 <- fit$beta_q975

    beta_resid     <- beta_mean - beta_true
    beta_ise       <- mean(beta_resid^2)
    beta_ise_rel   <- beta_ise / mean(beta_true^2)
    beta_bias_mean <- mean(beta_resid)
    beta_bias_max  <- max(abs(beta_resid))
    beta_cover_pw  <- mean((beta_true >= beta_q025) & (beta_true <= beta_q975))

    # Simultaneous coverage: is the WHOLE curve inside the max-statistic band?
    # Pointwise coverage is not a statement about the function.
    beta_cover_sim <- all(beta_true >= fit$beta_cma_lower &
                          beta_true <= fit$beta_cma_upper)
    beta_cma_qd    <- fit$beta_cma_qd

    ###############################################################
    ## scalar parameter metrics
    ###############################################################
    # gamma_true = DGM scale; gamma_implied = AFT target, and what bias/SE/
    # coverage are computed against.  They differ only for the Cox DGM.
    gamma_true    <- gamma_dgm
    gamma_implied <- if (is_cox) {
      aft_gamma_target(gamma_true, rho_dgm, LAMBDA_0)
    } else {
      gamma_true
    }
    gamma_est   <- colMeans(fit$gamma)
    gamma_bias  <- gamma_est - gamma_implied
    gamma_se    <- (gamma_est - gamma_implied)^2
    q_gamma     <- apply(fit$gamma, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    gamma_cover <- (gamma_implied >= q_gamma[1, ]) & (gamma_implied <= q_gamma[2, ])

    # tau2_implied = tau2_true/rho^2 for Cox, which given tau_dgm = rho*tau
    # evaluates to tau^2 -- the same target as the AFT arms.
    tau2_true    <- tau_dgm^2
    tau2_implied <- if (is_cox) tau2_true / rho_dgm^2 else tau2_true
    # Posterior MEDIAN, not mean: the tau2 posterior is right-skewed (skewness
    # 1.2-1.4 at J = 25), so the mean sits above the truth by construction.  At
    # tau = 1, J = 25 the median gives 1.119 against a mean of 1.188.
    tau2_est   <- median(fit$tau2)
    tau2_bias  <- tau2_est - tau2_implied
    tau2_se    <- (tau2_est - tau2_implied)^2
    tau2_cover <- (tau2_implied >= quantile(fit$tau2, 0.025)) &
      (tau2_implied <= quantile(fit$tau2, 0.975))

    # The target is the log-time RESIDUAL variance -- what the Gaussian fit
    # estimates -- not the DGM's scale parameter.  Both evaluate to sigma^2 by
    # the sigma_dgm construction, asserted outside the loop.
    sigma2_true  <- if (is_cox) pi^2 / (6 * rho_dgm^2) else
                    switch(as.character(family),
                           lognormal   = sigma_dgm^2,
                           loglogistic = sigma_dgm^2 * pi^2 / 3)
    sigma2_est   <- mean(fit$sigma2)
    sigma2_bias  <- sigma2_est - sigma2_true
    sigma2_se    <- (sigma2_est - sigma2_true)^2
    sigma2_cover <- (sigma2_true >= quantile(fit$sigma2, 0.025)) &
      (sigma2_true <= quantile(fit$sigma2, 0.975))

    ###############################################################
    ## out-of-sample predictive metrics -- MARGINAL (new clusters)
    ###############################################################
    pred_surv <- predict_gibbs_frailty(
      fit, X_new = test_data$X, Z_new = Z_test, cluster_id_new = NULL,
      type = "survival", times = tgrid, level = 0.95, quantiles = FALSE)

    # Negative posterior-mean linear predictor.  Rank-equivalent to -rowSums(S)
    # for an AFT model but cannot underflow: a saturated survival matrix
    # collapses -rowSums() to one repeated value, and UnoC scores only strict
    # lp_i > lp_j, so ties return C = 0 rather than 0.5.
    pred_link <- predict_gibbs_frailty(
      fit, X_new = test_data$X, Z_new = Z_test, cluster_id_new = NULL,
      type = "link", level = 0.95, quantiles = FALSE)
    risk_score <- -pred_link$mean

    Surv_train <- survival::Surv(data$Y, data$delta)
    c_index <- survAUC::UnoC(Surv.rsp = Surv_train, Surv.rsp.new = Surv_test,
                             lpnew = risk_score, time = tau_eval)
    ibs <- as.numeric(cal_IPCW_Brier(
      S_mat = pred_surv$mean, time_test = test_data$Y, event_test = test_data$delta,
      time_train = data$Y, event_train = data$delta, tgrid = tgrid,
      G_data = "test", eps_G = eps_G))

    ###############################################################
    ## out-of-sample predictive metrics -- CONDITIONAL (known clusters)
    ##
    ## New subjects in the SAME clusters as the training data, sharing its u_j
    ## via the DGM's `u` argument and held out of the fit.  Scored twice with the
    ## same fit and subjects: *_known uses u_hat_j, *_known_marg sets frailty to
    ## 0; the difference is the frailty's predictive contribution.  No leakage --
    ## these subjects' outcomes never entered the sampler.  Unlike the marginal
    ## test set this one cannot be fixed across iterations, so expect more noise.
    ###############################################################
    sim_known <- sim_dgm(N_known, n_cluster, u = sim_data$u)   # share the frailties
    known_data <- sim_known$data
    Z_known    <- model.matrix(~ Z1 + Z2, data = known_data)
    Surv_known <- survival::Surv(known_data$Y, known_data$delta)

    # predict_gibbs_frailty() silently falls back to frailty = 0 for any cluster
    # id it cannot match, which would relabel marginal predictions as
    # conditional.  Fail loudly instead.
    stopifnot(all(as.character(known_data$cluster_id) %in% fit$meta$cluster_levels))

    km_known    <- survival::survfit(survival::Surv(known_data$Y, 1 - known_data$delta) ~ 1)
    G_known     <- stats::stepfun(km_known$time, c(1, km_known$surv))
    tgrid_known <- tgrid[G_known(tgrid) >= eps_G]
    stopifnot(length(tgrid_known) >= 2L)
    tau_known   <- max(tgrid_known)

    eval_known <- function(cid) {
      lk <- predict_gibbs_frailty(fit, X_new = known_data$X, Z_new = Z_known,
              cluster_id_new = cid, type = "link", level = 0.95, quantiles = FALSE)
      sv <- predict_gibbs_frailty(fit, X_new = known_data$X, Z_new = Z_known,
              cluster_id_new = cid, type = "survival", times = tgrid_known,
              level = 0.95, quantiles = FALSE)
      c(c_index = survAUC::UnoC(Surv.rsp = Surv_train, Surv.rsp.new = Surv_known,
                                lpnew = -lk$mean, time = tau_known),
        ibs = as.numeric(cal_IPCW_Brier(
          S_mat = sv$mean, time_test = known_data$Y, event_test = known_data$delta,
          time_train = data$Y, event_train = data$delta, tgrid = tgrid_known,
          G_data = "test", eps_G = eps_G)))
    }
    m_cond <- eval_known(known_data$cluster_id)   # uses u_hat_j
    m_marg <- eval_known(NULL)                    # frailty forced to 0

    # How well is the frailty itself recovered?  fit$u columns follow
    # meta$cluster_levels, so index the true u through it rather than assuming
    # the orders agree.
    u_hat  <- colMeans(fit$u)
    u_true <- sim_data$u[as.integer(fit$meta$cluster_levels)]
    # For the Cox DGM u is on the log-hazard scale and u_AFT = -u_Cox/rho, so put
    # the truth on the AFT scale before correlating.
    if (is_cox) u_true <- -u_true / rho_dgm
    cor_u  <- suppressWarnings(cor(u_hat, u_true))

    ###############################################################
    ## COMPARATOR -- linear functional Cox with a shared frailty (fcoxFr)
    ##
    ## FPCA on X(s), then coxme with (1|cluster).  A frequentist, non-Bayesian
    ## alternative fitted to the SAME training data and scored on the SAME test
    ## sets, so lf_c_index / lf_ibs pair directly with c_index / ibs.
    ##
    ## Three things it does NOT give, which is why its columns are a subset:
    ##   * no interval estimates for beta(s), so no coverage
    ##   * beta_hat is on the LOG-HAZARD scale, so only its shape (correlation
    ##     with the truth, which is scale free) is comparable across families
    ##   * fpcr_predict() returns the FIXED-effect predictor only, so the
    ##     conditional score has to add ranef() by hand
    ###############################################################
    tic()
    fit_lf <- fpcr(time = data$Y, event = data$delta, group = data$cluster_id,
                   X = unclass(data$X),          # fda rejects an AsIs matrix
                   Z = Z[, -1, drop = FALSE], nb = 20, gp = s_grid)
    ts_lf   <- toc(quiet = TRUE)
    time_lf <- ts_lf$toc - ts_lf$tic

    # beta_hat is on the log-hazard scale; beta_AFT = -beta_Cox/rho, so its shape
    # matches -beta_true up to a positive factor.  Correlation is scale free and
    # therefore the one beta comparison valid for every family.
    lf_beta_cor <- suppressWarnings(cor(as.vector(fit_lf$bhat), -beta_true))
    lf_gamma_est <- unname(fit_lf$gammahat)
    lf_tau2_est  <- as.numeric(coxme::VarCorr(fit_lf$model)[[1]])

    lf_eta_train <- as.numeric(fit_lf$model$linear.predictor)
    lf_ranef     <- coxme::ranef(fit_lf$model)[[1]]

    lf_predict <- function(newdat, Znew, tg, cid = NULL) {
      pr <- fpcr_predict(fit_lf, time = newdat$Y, event = newdat$delta,
                         group = newdat$cluster_id, X = unclass(newdat$X),
                         Z = Znew[, -1, drop = FALSE], nb = 20, gp = s_grid)
      eta <- pr$predictions                       # fixed effects only
      if (!is.null(cid)) {
        j <- match(as.character(cid), names(lf_ranef))
        eta <- eta + ifelse(is.na(j), 0, lf_ranef[j])
      }
      S <- get_breslow_survival(lf_eta_train, data$Y, data$delta, eta, tg)
      c(c_index = survAUC::UnoC(Surv.rsp = Surv_train,
                                Surv.rsp.new = survival::Surv(newdat$Y, newdat$delta),
                                lpnew = eta, time = max(tg)),
        ibs = as.numeric(cal_IPCW_Brier(
          S_mat = S, time_test = newdat$Y, event_test = newdat$delta,
          time_train = data$Y, event_train = data$delta, tgrid = tg,
          G_data = "test", eps_G = eps_G)))
    }

    lf_marg      <- lf_predict(test_data, Z_test, tgrid)
    lf_cond_k    <- lf_predict(known_data, Z_known, tgrid_known, known_data$cluster_id)
    lf_marg_k    <- lf_predict(known_data, Z_known, tgrid_known)

    ###############################################################
    ## assemble
    ###############################################################
    beta_row <- list(mean = beta_mean, sd = beta_sd,
                     q025 = beta_q025, q975 = beta_q975)

    df_coef <- data.frame(
      beta_ise       = beta_ise,
      beta_ise_rel   = beta_ise_rel,   # ISE / mean(beta_true^2)
      beta_bias_mean = beta_bias_mean,
      beta_bias_max  = beta_bias_max,
      beta_cover_pw  = beta_cover_pw,
      beta_cover_sim = beta_cover_sim,
      beta_cma_qd    = beta_cma_qd,
      # frailty variance
      tau2_true    = tau2_true,
      tau2_implied = tau2_implied,
      tau2_est     = tau2_est,
      tau2_bias    = tau2_bias,
      tau2_se      = tau2_se,
      tau2_cover   = tau2_cover,
      # residual variance; = sigma^2 for both families by construction
      sigma2_true  = sigma2_true,
      sigma2_est   = sigma2_est,
      sigma2_bias  = sigma2_bias,
      sigma2_se    = sigma2_se,
      sigma2_cover = sigma2_cover,
      # predictive: discrimination (C-index) and calibration (IBS)
      c_index      = c_index,
      ibs          = ibs,
      # conditional, known clusters: *_known uses u_hat_j, *_known_marg frailty 0
      c_index_known      = m_cond[["c_index"]],
      ibs_known          = m_cond[["ibs"]],
      c_index_known_marg = m_marg[["c_index"]],
      ibs_known_marg     = m_marg[["ibs"]],
      cor_u              = cor_u,
      # a collapsed risk score drives Uno's C to 0 via ties; record enough to spot it
      risk_sd       = sd(risk_score),
      risk_n_unique = length(unique(risk_score)),
      # ---- comparator: linear functional Cox with frailty (fcoxFr) ----
      # Subset of the columns above -- see the block for why.  lf_* pairs with
      # the unprefixed column of the same name.
      lf_beta_cor        = lf_beta_cor,   # shape only; beta_hat is log-hazard scale
      lf_tau2_est        = lf_tau2_est,   # log-hazard scale frailty variance
      lf_c_index         = lf_marg[["c_index"]],
      lf_ibs             = lf_marg[["ibs"]],
      lf_c_index_known   = lf_cond_k[["c_index"]],
      lf_ibs_known       = lf_cond_k[["ibs"]],
      lf_c_index_known_marg = lf_marg_k[["c_index"]],
      lf_ibs_known_marg     = lf_marg_k[["ibs"]]
    )
    df_coef[paste0("lf_gamma_est_", seq_along(lf_gamma_est))] <- as.list(lf_gamma_est)

    idx <- seq_along(gamma_est) - 1L
    df_coef[paste0("gamma_true_",    idx)] <- as.list(gamma_true)
    df_coef[paste0("gamma_implied_", idx)] <- as.list(gamma_implied)
    df_coef[paste0("gamma_est_",     idx)] <- as.list(gamma_est)
    df_coef[paste0("gamma_bias_",  idx)] <- as.list(gamma_bias)
    df_coef[paste0("gamma_se_",    idx)] <- as.list(gamma_se)
    df_coef[paste0("gamma_cover_", idx)] <- as.list(gamma_cover)

    df_info <- data.frame(
      scenario    = scenario,
      iter        = iter,
      seed        = seeds[[iter]],
      family      = as.character(family),
      N_total     = N_total,
      n_cluster   = n_cluster,
      nS          = nS,
      beta_type   = as.character(beta_type),
      tau         = tau,          # AFT-scale target, common to all families
      tau_dgm     = tau_dgm,      # SD handed to the DGM (= rho * tau for Cox)
      sigma       = sigma,        # log-time residual SD target
      sigma_dgm   = sigma_dgm,    # scale parameter handed to the AFT DGM
      rho         = rho_dgm,      # Weibull shape used by the Cox DGM (NA otherwise)
      psi         = if (as.character(family) == "cox_tv") PSI_TV else 0,
      censor_rate = censor_rate,
      event_rate  = mean(data$delta),
      time        = time,        # Gibbs seconds
      time_lf     = time_lf      # comparator seconds
    )

    list(info = df_info, coef = df_coef, beta_row = beta_row)

  }, error = function(e) {
    warning(sprintf("Iteration %d skipped due to error:\n  %s", iter, e$message))
    NULL
  })

  if (!is.null(res)) {
    info_list[[iter]] <- res$info
    coef_list[[iter]] <- res$coef
    for (k in names(bmat)) bmat[[k]][iter, ] <- res$beta_row[[k]]
  }
} # end for loop

# Row i of every matrix is iteration i; use result$coef$iter to select the
# iterations that completed.  s_grid describes the columns of beta$*.
result <- list(
  info   = dplyr::bind_rows(Filter(Negate(is.null), info_list)),
  coef   = dplyr::bind_rows(Filter(Negate(is.null), coef_list)),
  s_grid = s_grid,
  beta_true = beta_true,
  beta   = bmat
)

###############################################################
## save result
###############################################################
Date = gsub("-", "", Sys.Date())
dir.create(file.path(here::here("outputs"), Date), showWarnings = FALSE)

filename = paste0(here::here("outputs", Date), "/", scenario, ".RDA")
save(result, file = filename)

###############################################################
## end sim
###############################################################
