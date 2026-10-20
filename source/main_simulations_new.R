####################################################################
# Weijia Qian
# Oct 11, 2026
#
# Simulate survival data from the EXTENDED HAZARDS bridge (simulate_EH), fit the
# functional frailty AFT model (gibbs_functional_frailty) and the linear
# functional Cox comparator (fcoxFr), and extract estimation + predictive metrics.
#
# This is the bridge-DGM counterpart of main_simulations.R and follows the same
# structure: design grid, cell-seeded test set shared across censoring levels,
# per-iteration seeds, known-cluster test set, comparator block, and the same
# result list saved to outputs/<Date>_EH/<scenario>.RDA.
#
# DGM:  h(t | a) = h0(t e^{kappa a}) e^{a},   a = -(z'gamma + int X beta + u_j),
#       h0 log-normal (SD of log T0 = sigma), u_j ~ N(0, tau^2)
#   kappa = 1   : exactly the functional frailty AFT model (model (1))
#   kappa = 0   : exactly a functional PH model with Gaussian frailty (fcoxFr's
#                 model); the baseline is not Weibull, so it has no AFT form
#   kappa = 0.5 : neither -- both working models are misspecified
#
# Targets.  Under misspecification there is no parameter value common to both
# models (the AFT beta is on the log-time scale, the Cox beta on the log-hazard
# scale, and for 0 < kappa < 1 neither scale is exact).  Two kinds of metric are
# therefore reported:
#   * TRUTH-BASED, SCALE-FREE -- comparable across the models at every kappa:
#       surv_iae / surv_rmse  error of the predicted survival curve against the
#                             TRUE curve, which the DGM knows exactly (no IPCW),
#                             PROBABILITY-weighted over the 5th-80th percentiles
#                             of the test set's true event-time distribution F:
#                               iae = mean_i (1/0.75) int_{0.05}^{0.80}
#                                       |S_hat_i(F^-1(p)) - S_i(F^-1(p))| dp,
#                             rmse analogous with squares.  Every percentile
#                             counts equally, whatever the time scale.
#       surv_iae_early/_late  the same over the 5th-40th percentiles (inside
#                             follow-up at every censoring level) and the
#                             40th-80th (beyond follow-up at heavy censoring, so
#                             it measures EXTRAPOLATION; fcoxFr's Breslow baseline
#                             is flat after the last training event)
#       beta_shape            || b_hat/||b_hat|| - b/||b|| ||^2 with each model's
#                             EXPECTED sign (+ AFT, - Cox): both models are
#                             single-index in int X beta, so the shape is their
#                             common estimand.  0 = perfect, 4 = sign-flipped
#       cor_u                 cor(predicted cluster effect, true u_j); for Cox the
#                             BLUPs estimate -u_j, so -BLUP is used
#   * EXACT -- against the true parameter, defined only where the truth lies in
#     the model's structural class: the AFT model at kappa = 1, Cox at kappa = 0.
#     NA elsewhere, deliberately (beta_ise, gamma_bias, tau2_bias, sigma2_bias,
#     coverage; lf_beta_ise, lf_gamma_bias, lf_tau2_bias).
# Secondary, as in main_simulations.R: Uno's IPCW C-index and the IPCW integrated
# Brier score, with ORACLE values from the true survival function as the ceiling.
#
# Time grids -- quantiles of the test set's TRUE event times, so common across
# censoring levels by construction:
#   tgrid_surv  p = 0.05, 0.075, ..., 0.80 (31 points): survival-curve error
#   tgrid       its points with p <= 0.40 (15 points): Uno's C at max(tgrid) and
#               the IBS, whose IPCW weights need G(t) > 0 (inside follow-up)
#
# Test sets and predictions
#   *  (no suffix)    NEW clusters, fixed per cell, frailty set to 0 (as in
#                     main_simulations.R); truth S(t | x, u = 0)
#   *_known           KNOWN clusters, redrawn per iteration, using u_hat_j;
#                     truth S(t | x, u_j)
#   *_known_marg      same subjects with the frailty forced to 0; truth
#                     S(t | x, u = 0)
#
# Per-iteration output (result$coef)
# ----------------------------------
#   beta_shape                              scale-free shape error (all kappa)
#   beta_{ise,ise_rel,bias_mean,bias_max}   exact; kappa = 1 only
#   beta_cover_pw / beta_cover_sim          exact; kappa = 1 only
#   beta_cma_qd                             simultaneous critical value
#   gamma_{true,est,bias,se,cover}_k        bias/se/cover exact; kappa = 1 only
#   tau2_{true,est,bias,se,cover}           bias/se/cover exact; kappa = 1 only
#   sigma2_{true,est,bias,se,cover}         bias/se/cover exact; kappa = 1 only
#   c_index, ibs                            new clusters, frailty = 0
#   surv_{iae,iae_early,iae_late,rmse}      new clusters, frailty = 0
#   {c_index,ibs,surv_*}_known              known clusters, u_hat
#   {c_index,ibs,surv_*}_known_marg         same subjects, frailty = 0
#   surv_ae_p10, ..., surv_ae_p80[_known[_marg]]
#                                           pointwise mean absolute error of the
#                                           predicted survival at the 10th, 20th,
#                                           ..., 80th percentiles of the true event
#                                           times (lf_ prefix for fcoxFr)
#   oracle_{c_index,ibs}[_known]            the true survival function's scores
#                                           (S(t | x, u = 0) for new clusters,
#                                           S(t | x, u_j) for known clusters)
#   cor_u                                   cor(u_hat_j, u_j)
#   risk_sd, risk_n_unique                  guards against a collapsed risk score
#   lf_*                                    fcoxFr comparator: lf_beta_shape,
#                                           lf_cor_u, the predictive metrics, and
#                                           exact lf_beta_ise / lf_gamma_bias /
#                                           lf_tau2_bias at kappa = 0 only
#
# result$info carries the design, event rate and sampler seconds; result$beta
# the fitted beta(s) on the evaluation grid; result$beta_true is beta(s) in
# model-(1) orientation (the exact AFT target at kappa = 1; the exact
# log-hazard target at kappa = 0 is -beta_true).
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
source(here("source", "simulate_EH.R"))
source(here("source", "fcoxFr_supp.R"))   # comparator: fpcr / fpcr_predict

###############################################################
## set simulation design elements
###############################################################
kappa = c(0, 0.5, 1)    # 0 = PH (Cox exact), 1 = AFT (AFT exact), 0.5 = neither
N_total = c(1000, 2000, 4000)
n_cluster = c(50)
nS = c(100)
beta_type = c("peak1")
tau = c(0.2, 0.5, 1)    # frailty SD on the log-time scale (exact target at both ends)
sigma = c(1)            # SD of the log-normal baseline's log T0
censor_rate = c(0.25, 0.5, 0.75)

N_iter = 500
N_test = 1000      # out-of-sample test set size (new subjects, NEW clusters)
N_known = 1000     # out-of-sample test set size (new subjects, KNOWN clusters);
                   # drawn per iteration sharing the training data's frailties

gamma_aft <- c(0.5, 0.3, -0.2)      # intercept, Z1, Z2 -- model-(1) orientation

# Fitting configuration: identical to main_simulations.R.
K_fit      <- 15
bs_fit     <- "cr"
n_iter_fit <- 20000
n_burn_fit <- 10000

params = expand.grid(kappa = kappa,
                     N_total = N_total,
                     n_cluster = n_cluster,
                     nS = nS,
                     beta_type = beta_type,
                     tau = tau,
                     sigma = sigma,
                     censor_rate = censor_rate)

## define number of simulations and parameter scenarios
if(doLocal) {
  scenario = 2
  N_iter = 1
}else{
  # defined from batch script params
  scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
}

# Validate here rather than letting an out-of-range or NA index propagate.
# nrow(params) changes whenever the design grid changes.
if (length(scenario) != 1L || is.na(scenario) ||
    scenario < 1L || scenario > nrow(params)) {
  stop(sprintf(
    "scenario must be a single integer in 1..%d, got: %s\n  (doLocal = %s, commandArgs = %s)",
    nrow(params), paste(deparse(scenario), collapse = ""), doLocal,
    paste(deparse(commandArgs(trailingOnly = TRUE)), collapse = "")))
}

kappa       = params$kappa[scenario]
N_total     = params$N_total[scenario]
n_cluster   = params$n_cluster[scenario]
nS          = params$nS[scenario]
beta_type   = params$beta_type[scenario]
tau         = params$tau[scenario]
sigma       = params$sigma[scenario]
censor_rate = params$censor_rate[scenario]

# Exact targets exist only where the truth lies in the model's structural class.
aft_exact <- isTRUE(all.equal(kappa, 1))
cox_exact <- isTRUE(all.equal(kappa, 0))

# One DGM call.  `u` and N are the only things that vary between the training,
# marginal-test and known-cluster draws.
sim_dgm <- function(N, J, u = NULL) {
  simulate_EH(N_total = N, n_cluster = J, u = u, kappa = kappa,
              baseline = "lognormal", nS = nS,
              beta_type = as.character(beta_type), gamma = gamma_aft,
              sigma = sigma, tau = tau, censor_rate = censor_rate)
}

###############################################################
## metric helpers
###############################################################
# Error of predicted against TRUE survival (rows = subjects, cols = the time
# points t_k = F^-1(p_k)), PROBABILITY-weighted: trapezoidal in p and normalised
# by the width of the p-window, so every percentile of the event-time
# distribution counts equally.  Weighting in t instead would hand most of the
# weight to the long right tail (78-84% to the 40th-80th percentiles).
surv_err <- function(S_hat, S_true, p) {
  w <- .eh_trapz_weights(p) / (max(p) - min(p))
  d <- S_hat - S_true
  c(iae = mean(abs(d) %*% w), rmse = sqrt(mean((d^2) %*% w)))
}

# Whole window (5th-80th percentile) and its early (5th-40th) and late
# (40th-80th) parts.  p_surv, i_early and i_late are set with the test set below.
surv_metrics <- function(S_hat, S_true) {
  part <- function(i) surv_err(S_hat[, i, drop = FALSE], S_true[, i, drop = FALSE],
                               p_surv[i])[["iae"]]
  all <- surv_err(S_hat, S_true, p_surv)
  # pointwise mean absolute error at the 10th, 20th, ..., 80th percentiles
  ae  <- colMeans(abs(S_hat[, i_pct, drop = FALSE] - S_true[, i_pct, drop = FALSE]))
  c(surv_iae = all[["iae"]], surv_iae_early = part(i_early),
    surv_iae_late = part(i_late), surv_rmse = all[["rmse"]],
    setNames(ae, paste0("surv_ae_p", round(100 * p_surv[i_pct]))))
}

# Scale-free shape error with a FIXED expected sign (0 = same shape, 4 = flipped).
shape_err <- function(b_hat, b_true, sgn) {
  nrm <- function(b) b / sqrt(mean(b^2))
  if (!all(is.finite(b_hat)) || sqrt(mean(b_hat^2)) < 1e-12) return(NA_real_)
  mean((nrm(sgn * b_hat) - nrm(b_true))^2)
}

###############################################################
## run simulations
###############################################################

coef_list <- vector("list", length = N_iter)
info_list <- vector("list", length = N_iter)

###############################################################
## Marginal test set: one draw shared by every iteration.
## New clusters -> frailty = 0 at prediction time.
##
## Seeded on the CELL, not the scenario.  The cell is every design factor except
## censor_rate, so the three censoring levels of a cell share one seed.  Because
## simulate_EH consumes censor_rate only AFTER drawing T_true, that makes T_true
## byte-identical across the three, and the censoring comparison is PAIRED.
###############################################################
cell_cols <- setdiff(names(params), "censor_rate")
cell_seed <- match(do.call(paste, c(params[scenario, cell_cols, drop = FALSE], sep = "|")),
                   unique(do.call(paste, c(params[, cell_cols, drop = FALSE], sep = "|"))))
stopifnot(length(cell_seed) == 1L, !is.na(cell_seed))
set.seed(cell_seed)

# Deliberately NOT a mirror of the training design.  Predictions here set the
# frailty to 0, so the test clusters exist only to inject between-cluster
# variance; since this set is drawn ONCE and reused by every iteration, few
# clusters would mean a small sample of u_j whose idiosyncrasy biases the whole
# scenario and never averages out.
n_cluster_test <- max(50L, N_test %/% 5L)

sim_test  <- sim_dgm(N_test, n_cluster_test)
test_data <- sim_test$data
Z_test    <- model.matrix(~ Z1 + Z2, data = test_data)
s_grid    <- sim_test$s_grid
cat("Test set generated: N =", nrow(test_data),
    "| event rate =", round(mean(test_data$delta), 3), "\n")

# Two nested grids, both quantiles of the TRUE event-time distribution of the
# test set, and so common across censoring levels by construction.
#   tgrid_surv  p = 0.05, 0.075, ..., 0.80: survival-curve error.  Scored against
#               the true curve, it needs no IPCW and may run past follow-up;
#               there it measures extrapolation.
#   tgrid       the points with p <= HORIZON_Q: Uno's C and the IBS.  Their IPCW
#               weights need G(t) > 0, and uniform censoring has bounded support,
#               so this grid stops inside follow-up even at 75% censoring.
# tgrid is a subset of tgrid_surv so that ONE prediction serves both: every
# extra time point costs the AFT model a pass over all posterior draws.
HORIZON_Q <- 0.40
SURV_Q    <- 0.80
p_surv     <- seq(0.05, SURV_Q, by = 0.025)
tgrid_surv <- unname(quantile(test_data$T_true, probs = p_surv))
stopifnot(!anyDuplicated(tgrid_surv))
i_early <- which(p_surv <= HORIZON_Q + 1e-9)    # 5th-40th percentile
i_late  <- which(p_surv >= HORIZON_Q - 1e-9)    # 40th-80th percentile
i_ipcw  <- i_early                              # columns of tgrid_surv in tgrid
i_pct   <- match(round(seq(0.10, SURV_Q, by = 0.10), 3), round(p_surv, 3))
stopifnot(!anyNA(i_pct))                        # 10th, 20th, ..., 80th percentiles

eps_G      <- 0.05
km_cens_te <- survival::survfit(survival::Surv(test_data$Y, 1 - test_data$delta) ~ 1)
G_test     <- stats::stepfun(km_cens_te$time, c(1, km_cens_te$surv))
if (any(G_test(tgrid_surv[i_ipcw]) < eps_G)) {
  warning(sprintf(paste("Scenario %d: G(t) < %.2f at %d of %d IPCW grid points;",
                        "the horizon is NO LONGER common across censoring levels.",
                        "Lower HORIZON_Q (currently %.2f)."),
                  scenario, eps_G, sum(G_test(tgrid_surv[i_ipcw]) < eps_G),
                  length(i_ipcw), HORIZON_Q))
  i_ipcw <- i_ipcw[G_test(tgrid_surv[i_ipcw]) >= eps_G]
}
stopifnot(length(i_ipcw) >= 2L)
tgrid <- tgrid_surv[i_ipcw]

tau_eval <- max(tgrid)
cat("Evaluation horizon: tgrid has", length(tgrid), "points, tau =",
    round(tau_eval, 4), "| G(tau) =", round(G_test(tau_eval), 4),
    "| survival grid to", round(max(tgrid_surv), 4), "\n")

Surv_test <- survival::Surv(test_data$Y, test_data$delta)

# TRUE survival of the test subjects at frailty 0 -- the target of a
# frailty = 0 prediction.
S_test_true <- sim_test$surv_true(tgrid_surv, test_data$lp, 0)

# beta(s) in model-(1) orientation, on the simulation grid.
beta_true <- sim_test$beta
n_point   <- length(s_grid)

mat_grid <- function() matrix(NA_real_, N_iter, n_point)
bmat     <- list(mean = mat_grid(), sd = mat_grid(),
                 q025 = mat_grid(), q975 = mat_grid())

# Per-iteration seeds
# Drawn from the cell stream seeded above, so the three censoring levels of a
# cell also share their per-iteration training seeds: same subjects, same true
# event times, only the censoring differs.
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
      K = K_fit, basis_type = bs_fit,
      lambda_init = 1, A_lambda = 1, B_lambda = 1,
      var_gamma = 100,
      A_tau2 = 0.01, B_tau2 = 0.01,
      A_sigma2 = 0.01, B_sigma2 = 0.01, 
      n_iter = n_iter_fit, n_burn = n_burn_fit, n_thin = 1, verbose = TRUE
    )
    time_stamp <- toc(quiet = TRUE)
    time <- time_stamp$toc - time_stamp$tic

    ###############################################################
    ## beta(s) estimation metrics
    ###############################################################
    beta_mean <- fit$beta_mean
    beta_sd   <- fit$beta_cma_sd
    beta_q025 <- fit$beta_q025
    beta_q975 <- fit$beta_q975

    # Scale-free shape error: defined at every kappa.
    beta_shape <- shape_err(beta_mean, beta_true, +1)

    # Exact metrics: the AFT coefficient function is the true one only at kappa = 1.
    beta_resid     <- beta_mean - beta_true
    beta_ise       <- if (aft_exact) mean(beta_resid^2) else NA_real_
    beta_ise_rel   <- if (aft_exact) beta_ise / mean(beta_true^2) else NA_real_
    beta_bias_mean <- if (aft_exact) mean(beta_resid) else NA_real_
    beta_bias_max  <- if (aft_exact) max(abs(beta_resid)) else NA_real_
    beta_cover_pw  <- if (aft_exact) mean((beta_true >= beta_q025) & (beta_true <= beta_q975)) else NA_real_
    beta_cover_sim <- if (aft_exact) all(beta_true >= fit$beta_cma_lower &
                                         beta_true <= fit$beta_cma_upper) else NA
    beta_cma_qd    <- fit$beta_cma_qd

    ###############################################################
    ## scalar parameter metrics (exact at kappa = 1 only)
    ###############################################################
    gamma_true  <- gamma_aft
    gamma_est   <- colMeans(fit$gamma)
    q_gamma     <- apply(fit$gamma, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    gamma_bias  <- if (aft_exact) gamma_est - gamma_true else rep(NA_real_, length(gamma_est))
    gamma_se    <- if (aft_exact) (gamma_est - gamma_true)^2 else rep(NA_real_, length(gamma_est))
    gamma_cover <- if (aft_exact) (gamma_true >= q_gamma[1, ]) & (gamma_true <= q_gamma[2, ]) else
                   rep(NA, length(gamma_est))

    # Posterior MEDIAN: the tau2 posterior is right-skewed.
    tau2_true  <- tau^2
    tau2_est   <- median(fit$tau2)
    tau2_bias  <- if (aft_exact) tau2_est - tau2_true else NA_real_
    tau2_se    <- if (aft_exact) (tau2_est - tau2_true)^2 else NA_real_
    tau2_cover <- if (aft_exact) (tau2_true >= quantile(fit$tau2, 0.025)) &
                                 (tau2_true <= quantile(fit$tau2, 0.975)) else NA

    sigma2_true  <- sigma^2                   # log-time residual variance at kappa = 1
    sigma2_est   <- mean(fit$sigma2)
    sigma2_bias  <- if (aft_exact) sigma2_est - sigma2_true else NA_real_
    sigma2_se    <- if (aft_exact) (sigma2_est - sigma2_true)^2 else NA_real_
    sigma2_cover <- if (aft_exact) (sigma2_true >= quantile(fit$sigma2, 0.025)) &
                                   (sigma2_true <= quantile(fit$sigma2, 0.975)) else NA

    ###############################################################
    ## out-of-sample predictive metrics -- MARGINAL (new clusters)
    ###############################################################
    pred_surv <- predict_gibbs_frailty(
      fit, X_new = test_data$X, Z_new = Z_test, cluster_id_new = NULL,
      type = "survival", times = tgrid_surv, level = 0.95, quantiles = FALSE)$mean

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
      S_mat = pred_surv[, i_ipcw, drop = FALSE], time_test = test_data$Y,
      event_test = test_data$delta, time_train = data$Y, event_train = data$delta,
      tgrid = tgrid, G_data = "test", eps_G = eps_G))
    e_test <- surv_metrics(pred_surv, S_test_true)

    # Oracle: the true survival function scored the same way.
    oracle_c_index <- survAUC::UnoC(Surv.rsp = Surv_train, Surv.rsp.new = Surv_test,
                                    lpnew = -rowSums(S_test_true[, i_ipcw, drop = FALSE]),
                                    time = tau_eval)
    oracle_ibs <- as.numeric(cal_IPCW_Brier(
      S_mat = S_test_true[, i_ipcw, drop = FALSE], time_test = test_data$Y,
      event_test = test_data$delta,
      time_train = data$Y, event_train = data$delta, tgrid = tgrid,
      G_data = "test", eps_G = eps_G))

    ###############################################################
    ## out-of-sample predictive metrics -- CONDITIONAL (known clusters)
    ##
    ## New subjects in the SAME clusters as the training data, sharing its u_j
    ## via the DGM's `u` argument and held out of the fit.  Scored twice with the
    ## same fit and subjects: *_known uses u_hat_j against S(t | x, u_j);
    ## *_known_marg sets the frailty to 0 against S(t | x, u = 0).
    ###############################################################
    sim_known <- sim_dgm(N_known, n_cluster, u = sim_data$u)   # share the frailties
    known_data <- sim_known$data
    Z_known    <- model.matrix(~ Z1 + Z2, data = known_data)
    Surv_known <- survival::Surv(known_data$Y, known_data$delta)

    # predict_gibbs_frailty() silently falls back to frailty = 0 for any cluster
    # id it cannot match, which would relabel marginal predictions as
    # conditional.  Fail loudly instead.
    stopifnot(all(as.character(known_data$cluster_id) %in% fit$meta$cluster_levels))

    # Hold the horizon fixed; warn instead of silently shortening.
    km_known    <- survival::survfit(survival::Surv(known_data$Y, 1 - known_data$delta) ~ 1)
    G_known     <- stats::stepfun(km_known$time, c(1, km_known$surv))
    tgrid_known <- tgrid
    tau_known   <- tau_eval
    if (any(G_known(tgrid_known) < eps_G))
      warning(sprintf("Iter %d: known-set G(t) < %.2f at %d of %d grid points.",
                      iter, eps_G, sum(G_known(tgrid_known) < eps_G), length(tgrid_known)))

    # True survival of the known-cluster subjects on the survival grid (the
    # IPCW metrics use its first columns, i_ipcw).
    S_known_true      <- sim_known$surv_true(tgrid_surv, known_data$lp,
                                             sim_data$u[known_data$cluster_id])
    S_known_marg_true <- sim_known$surv_true(tgrid_surv, known_data$lp, 0)

    eval_known <- function(cid, S_truth) {
      lk <- predict_gibbs_frailty(fit, X_new = known_data$X, Z_new = Z_known,
              cluster_id_new = cid, type = "link", level = 0.95, quantiles = FALSE)
      sv <- predict_gibbs_frailty(fit, X_new = known_data$X, Z_new = Z_known,
              cluster_id_new = cid, type = "survival", times = tgrid_surv,
              level = 0.95, quantiles = FALSE)
      c(c_index = survAUC::UnoC(Surv.rsp = Surv_train, Surv.rsp.new = Surv_known,
                                lpnew = -lk$mean, time = tau_known),
        ibs = as.numeric(cal_IPCW_Brier(
          S_mat = sv$mean[, i_ipcw, drop = FALSE], time_test = known_data$Y,
          event_test = known_data$delta, time_train = data$Y,
          event_train = data$delta, tgrid = tgrid_known,
          G_data = "test", eps_G = eps_G)),
        surv_metrics(sv$mean, S_truth))
    }
    m_cond <- eval_known(known_data$cluster_id, S_known_true)    # uses u_hat_j
    m_marg <- eval_known(NULL, S_known_marg_true)                # frailty forced to 0

    oracle_c_index_known <- survAUC::UnoC(Surv.rsp = Surv_train, Surv.rsp.new = Surv_known,
                                          lpnew = -rowSums(S_known_true[, i_ipcw, drop = FALSE]),
                                          time = tau_known)
    oracle_ibs_known <- as.numeric(cal_IPCW_Brier(
      S_mat = S_known_true[, i_ipcw, drop = FALSE], time_test = known_data$Y,
      event_test = known_data$delta,
      time_train = data$Y, event_train = data$delta, tgrid = tgrid_known,
      G_data = "test", eps_G = eps_G))

    # How well is the frailty itself recovered?  fit$u columns follow
    # meta$cluster_levels, so index the true u through it.
    u_hat  <- colMeans(fit$u)
    u_true <- sim_data$u[as.integer(fit$meta$cluster_levels)]
    cor_u  <- suppressWarnings(cor(u_hat, u_true))

    ###############################################################
    ## COMPARATOR -- linear functional Cox with a shared frailty (fcoxFr)
    ##
    ## FPCA on X(s), then coxme with (1|cluster), fitted to the SAME training data
    ## and scored on the SAME test sets.  Its coefficients are on the log-hazard
    ## scale, so the exact targets (kappa = 0 only) are -beta(s), -gamma[-1] and
    ## tau^2; the shape and predictive metrics need no target and are defined at
    ## every kappa.
    ###############################################################
    tic()
    fit_lf <- fpcr(time = data$Y, event = data$delta, group = data$cluster_id,
                   X = unclass(data$X),          # fda rejects an AsIs matrix
                   Z = Z[, -1, drop = FALSE], nb = 20, gp = s_grid)
    ts_lf   <- toc(quiet = TRUE)
    time_lf <- ts_lf$toc - ts_lf$tic

    lf_bhat       <- as.vector(fit_lf$bhat)
    lf_beta_shape <- shape_err(lf_bhat, beta_true, -1)
    lf_beta_ise   <- if (cox_exact) mean((lf_bhat + beta_true)^2) else NA_real_

    lf_gamma_est  <- unname(fit_lf$gammahat)
    stopifnot(length(lf_gamma_est) == ncol(Z) - 1L)
    lf_gamma_bias <- if (cox_exact) lf_gamma_est + gamma_aft[-1] else
                     rep(NA_real_, length(lf_gamma_est))

    lf_tau2_est  <- as.numeric(coxme::VarCorr(fit_lf$model)[[1]])
    lf_tau2_bias <- if (cox_exact) lf_tau2_est - tau2_true else NA_real_

    # coxme's linear.predictor includes the random effects, so the Breslow
    # baseline below is the CONDITIONAL baseline.  It is a step function that
    # stays flat after the last training event, so beyond follow-up fcoxFr's
    # predicted survival stops falling -- the _late error measures that.
    lf_eta_train <- as.numeric(fit_lf$model$linear.predictor)
    lf_ranef     <- coxme::ranef(fit_lf$model)[[1]]
    lf_cor_u     <- suppressWarnings(cor(-lf_ranef[as.character(seq_len(n_cluster))], sim_data$u))
    lf_H0        <- -log(get_breslow_survival(lf_eta_train, data$Y, data$delta, 0, tgrid_surv))[1, ]

    lf_predict <- function(newdat, Znew, S_truth, cid = NULL) {
      pr <- fpcr_predict(fit_lf, time = newdat$Y, event = newdat$delta,
                         group = newdat$cluster_id, X = unclass(newdat$X),
                         Z = Znew[, -1, drop = FALSE], nb = 20, gp = s_grid)
      eta <- pr$predictions                       # fixed effects only
      if (!is.null(cid)) {
        j <- match(as.character(cid), names(lf_ranef))
        eta <- eta + ifelse(is.na(j), 0, lf_ranef[j])
      }
      S <- exp(-outer(exp(eta), lf_H0))           # = get_breslow_survival(lf_eta_train, ..., eta, tgrid_surv)
      c(c_index = survAUC::UnoC(Surv.rsp = Surv_train,
                                Surv.rsp.new = survival::Surv(newdat$Y, newdat$delta),
                                lpnew = eta, time = tau_eval),
        ibs = as.numeric(cal_IPCW_Brier(
          S_mat = S[, i_ipcw, drop = FALSE], time_test = newdat$Y,
          event_test = newdat$delta, time_train = data$Y,
          event_train = data$delta, tgrid = tgrid,
          G_data = "test", eps_G = eps_G)),
        surv_metrics(S, S_truth))
    }

    lf_marg   <- lf_predict(test_data, Z_test, S_test_true)
    lf_cond_k <- lf_predict(known_data, Z_known, S_known_true,
                            cid = known_data$cluster_id)
    lf_marg_k <- lf_predict(known_data, Z_known, S_known_marg_true)

    ###############################################################
    ## assemble
    ###############################################################
    beta_row <- list(mean = beta_mean, sd = beta_sd,
                     q025 = beta_q025, q975 = beta_q975)

    df_coef <- data.frame(
      beta_shape     = beta_shape,     # scale-free, all kappa
      beta_ise       = beta_ise,       # exact, kappa = 1 only
      beta_ise_rel   = beta_ise_rel,
      beta_bias_mean = beta_bias_mean,
      beta_bias_max  = beta_bias_max,
      beta_cover_pw  = beta_cover_pw,
      beta_cover_sim = beta_cover_sim,
      beta_cma_qd    = beta_cma_qd,
      # frailty variance (bias/se/cover exact, kappa = 1 only)
      tau2_true    = tau2_true,
      tau2_est     = tau2_est,
      tau2_bias    = tau2_bias,
      tau2_se      = tau2_se,
      tau2_cover   = tau2_cover,
      # residual variance (bias/se/cover exact, kappa = 1 only)
      sigma2_true  = sigma2_true,
      sigma2_est   = sigma2_est,
      sigma2_bias  = sigma2_bias,
      sigma2_se    = sigma2_se,
      sigma2_cover = sigma2_cover,
      # predictive, NEW clusters, frailty = 0
      c_index      = c_index,
      ibs          = ibs,
      surv_iae       = e_test[["surv_iae"]],         # 5th-80th percentile
      surv_iae_early = e_test[["surv_iae_early"]],   # 5th-40th
      surv_iae_late  = e_test[["surv_iae_late"]],    # 40th-80th
      surv_rmse      = e_test[["surv_rmse"]],
      oracle_c_index = oracle_c_index,
      oracle_ibs     = oracle_ibs,
      # predictive, KNOWN clusters: *_known uses u_hat_j, *_known_marg frailty 0
      c_index_known       = m_cond[["c_index"]],
      ibs_known           = m_cond[["ibs"]],
      surv_iae_known       = m_cond[["surv_iae"]],
      surv_iae_early_known = m_cond[["surv_iae_early"]],
      surv_iae_late_known  = m_cond[["surv_iae_late"]],
      surv_rmse_known      = m_cond[["surv_rmse"]],
      c_index_known_marg  = m_marg[["c_index"]],
      ibs_known_marg      = m_marg[["ibs"]],
      surv_iae_known_marg       = m_marg[["surv_iae"]],
      surv_iae_early_known_marg = m_marg[["surv_iae_early"]],
      surv_iae_late_known_marg  = m_marg[["surv_iae_late"]],
      surv_rmse_known_marg      = m_marg[["surv_rmse"]],
      oracle_c_index_known = oracle_c_index_known,
      oracle_ibs_known     = oracle_ibs_known,
      cor_u              = cor_u,
      # a collapsed risk score drives Uno's C to 0 via ties; record enough to spot it
      risk_sd       = sd(risk_score),
      risk_n_unique = length(unique(risk_score)),
      # ---- comparator: linear functional Cox with frailty (fcoxFr) ----
      lf_beta_shape  = lf_beta_shape,  # scale-free, all kappa
      lf_beta_ise    = lf_beta_ise,    # exact (log-hazard scale, -beta), kappa = 0 only
      lf_tau2_est    = lf_tau2_est,
      lf_tau2_bias   = lf_tau2_bias,   # exact, kappa = 0 only
      lf_cor_u       = lf_cor_u,
      lf_c_index     = lf_marg[["c_index"]],
      lf_ibs         = lf_marg[["ibs"]],
      lf_surv_iae       = lf_marg[["surv_iae"]],
      lf_surv_iae_early = lf_marg[["surv_iae_early"]],
      lf_surv_iae_late  = lf_marg[["surv_iae_late"]],
      lf_surv_rmse      = lf_marg[["surv_rmse"]],
      lf_c_index_known       = lf_cond_k[["c_index"]],
      lf_ibs_known           = lf_cond_k[["ibs"]],
      lf_surv_iae_known       = lf_cond_k[["surv_iae"]],
      lf_surv_iae_early_known = lf_cond_k[["surv_iae_early"]],
      lf_surv_iae_late_known  = lf_cond_k[["surv_iae_late"]],
      lf_surv_rmse_known      = lf_cond_k[["surv_rmse"]],
      lf_c_index_known_marg  = lf_marg_k[["c_index"]],
      lf_ibs_known_marg      = lf_marg_k[["ibs"]],
      lf_surv_iae_known_marg       = lf_marg_k[["surv_iae"]],
      lf_surv_iae_early_known_marg = lf_marg_k[["surv_iae_early"]],
      lf_surv_iae_late_known_marg  = lf_marg_k[["surv_iae_late"]],
      lf_surv_rmse_known_marg      = lf_marg_k[["surv_rmse"]]
    )
    # Pointwise errors at the 10th, ..., 80th percentiles, every model and setting.
    add_ae <- function(df, v, prefix = "", suffix = "") {
      a <- v[grep("^surv_ae_p", names(v))]
      df[paste0(prefix, names(a), suffix)] <- as.list(unname(a))
      df
    }
    df_coef <- add_ae(df_coef, e_test)
    df_coef <- add_ae(df_coef, m_cond,    suffix = "_known")
    df_coef <- add_ae(df_coef, m_marg,    suffix = "_known_marg")
    df_coef <- add_ae(df_coef, lf_marg,   prefix = "lf_")
    df_coef <- add_ae(df_coef, lf_cond_k, prefix = "lf_", suffix = "_known")
    df_coef <- add_ae(df_coef, lf_marg_k, prefix = "lf_", suffix = "_known_marg")

    df_coef[paste0("lf_gamma_est_",  seq_along(lf_gamma_est))]  <- as.list(lf_gamma_est)
    df_coef[paste0("lf_gamma_bias_", seq_along(lf_gamma_bias))] <- as.list(lf_gamma_bias)

    idx <- seq_along(gamma_est) - 1L
    df_coef[paste0("gamma_true_",  idx)] <- as.list(gamma_true)
    df_coef[paste0("gamma_est_",   idx)] <- as.list(gamma_est)
    df_coef[paste0("gamma_bias_",  idx)] <- as.list(gamma_bias)
    df_coef[paste0("gamma_se_",    idx)] <- as.list(gamma_se)
    df_coef[paste0("gamma_cover_", idx)] <- as.list(gamma_cover)

    df_info <- data.frame(
      scenario    = scenario,
      iter        = iter,
      seed        = seeds[[iter]],
      kappa       = kappa,        # 0 = PH, 1 = AFT, 0.5 = neither
      aft_exact   = aft_exact,
      cox_exact   = cox_exact,
      N_total     = N_total,
      n_cluster   = n_cluster,
      nS          = nS,
      beta_type   = as.character(beta_type),
      tau         = tau,
      sigma       = sigma,
      censor_rate = censor_rate,
      tau_eval    = tau_eval,      # IPCW horizon (40th pct); common across censoring
      n_tgrid     = length(tgrid),
      tau_surv    = max(tgrid_surv),  # end of the survival-error window (80th pct)
      n_tgrid_surv = length(tgrid_surv),
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
## outputs/<Date>_EH/ (as main_simulations_interaction.R uses <Date>_interaction):
## simulation_results.Rmd auto-loads the newest outputs/YYYYMMDD folder, so this
## study must not write there.
###############################################################
Date = gsub("-", "", Sys.Date())
dir.create(file.path(here::here("outputs"), paste0(Date, "_EH")), showWarnings = FALSE, recursive = TRUE)

filename = paste0(here::here("outputs", paste0(Date, "_EH")), "/", scenario, ".RDA")
save(result, file = filename)

###############################################################
## end sim
###############################################################
