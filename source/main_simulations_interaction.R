####################################################################
# Weijia Qian
# Jan 10, 2026
#
# Simulate survival data under a functional frailty AFT or Cox model with a
# bivariate surface beta(s, age), fit gibbs_frailty_interaction(), and extract
# estimation + predictive metrics.  Full derivations and the rationale for the
# scale corrections are in analysis/simulation_summary_interaction.md.
#
# Everything is scored on the AFT (log-time) scale.  For the Cox DGM the target
# is the implied AFT quantity, not the DGM parameter: *_implied columns hold the
# target, *_true the DGM-scale value.
#
# Per-iteration output (result$coef), and ni_* for the no-interaction comparator
# ----------------------------------------------------------------------------
#   beta_{ise,ise_rel,bias_mean,bias_max}   surface error; ise_rel = ise/mean(beta^2)
#   beta_cover_pw / beta_cover_sim          pointwise / simultaneous 95% coverage
#   beta_cma_qd                             simultaneous critical value
#   gamma_{true,implied,est,bias,se,cover}_k
#   tau2_{true,implied,est,bias,se,cover}
#   sigma2_{true,est,bias,se,cover}         target is the log-time residual variance
#   c_index, ibs                            marginal: new clusters, frailty = 0
#   {c_index,ibs}_known                     conditional: known clusters, frailty = u_hat
#   {c_index,ibs}_known_marg                same subjects, frailty = 0
#   cor_u                                   cor(u_hat_j, u_j) on the AFT scale
#   risk_sd, risk_n_unique                  guards against a collapsed risk score
#
# result$info carries the design, event rate and sampler seconds (time, time_ni);
# result$surface the fitted surface on the evaluation grid.
####################################################################

suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(splines))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survAUC))    # for UnoC (IPCW C-index)
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(tidyverse))

wd = getwd()

if(substring(wd, 2, 6) == "Users"){
  doLocal = TRUE
}else{
  doLocal = FALSE
}

###############################################################
## define or source functions used in code below
###############################################################
source(here("source", "gibbs.R"))
source(here("source", "gibbs_frailty_interaction.R"))
source(here("source", "predict_gibbs_interaction.R"))
# source(here("source", "predict_gibbs.R"))       # comparator predictions
source(here("source", "cal_Brier.R"))
source(here("source", "simulate_AFT_interaction.R"))
source(here("source", "simulate_Cox_interaction.R"))

###############################################################
## set simulation design elements
###############################################################
family = c("lognormal", "loglogistic", "cox")
N_total = c(1000, 2000, 4000)
n_cluster = c(25, 50, 75)
nS = c(100)
beta_type = c("additive")
# beta_type = c("additive", "bilinear", "ridge")
tau = c(0.2, 0.5, 1)  # AFT-scale (log-time) frailty SD target, common to all
                      # families; the Cox arm is simulated at rho * tau
sigma = c(1)       # AFT-scale (log-time) RESIDUAL SD target, common to all
                   # families.  Each DGM takes a different scale parameter to
                   # hit it; for the Cox arm it also fixes rho.  See below.
censor_rate = c(0.25, 0.5, 0.75)
N_iter = 500
N_test = 1000      # out-of-sample test set size (new subjects, NEW clusters)
N_known = 1000     # out-of-sample test set size (new subjects, KNOWN clusters);
                   # drawn per iteration sharing the training data's frailties

# Fitting configuration, shared by both models so that the only difference
# between them is whether the coefficient varies with age.
K_s_fit   <- 15
K_age_fit <- 10
bs_s_fit  <- "cr"      # "cc" for real NHANES data, "cr" for simulation
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
  scenario = 100
  N_iter = 1
}else{
  # defined from batch script params
  scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
}

family      = params$family[scenario]
N_total     = params$N_total[scenario]
n_cluster   = params$n_cluster[scenario]
nS          = params$nS[scenario]
beta_type   = params$beta_type[scenario]
tau         = params$tau[scenario]
sigma       = params$sigma[scenario]
censor_rate = params$censor_rate[scenario]

# sigma is the shared log-time RESIDUAL SD target, so all families are compared
# at matched signal-to-noise and only the error SHAPE differs.  Each DGM needs a
# different scale parameter to hit it, since sigma_err = sigma_dgm * sd(epsilon).
# For the Cox arm the log-time error is Gumbel_min(0,1)/rho, so rho is not free.
sigma_dgm <- switch(as.character(family),
                    lognormal   = sigma,
                    loglogistic = sigma * sqrt(3) / pi,
                    cox         = NA_real_)

rho_dgm <- if (as.character(family) == "cox") pi / (sqrt(6) * sigma) else NA_real_

# gamma is specified on the AFT scale.  gamma_AFT = -gamma_Cox/rho, so the Cox
# arm is handed rho*gamma to give its Z slopes the same AFT-scale magnitude.
# Only the slopes are matched; the intercept is a pure time-scale offset that
# the censoring calibration and the Brier t-grid both absorb.
gamma_aft <- c(0.5, 0.3, -0.2)      # intercept, Z1, Z2 -- AFT (log-time) scale
gamma_dgm <- if (as.character(family) == "cox") rho_dgm * gamma_aft else gamma_aft

# Guard: the derivation above must leave every family with the SAME log-time
# residual variance, sigma^2.  If this fires, sigma_dgm / rho_dgm and the
# sigma2_true switch inside the loop have drifted out of sync.
stopifnot(isTRUE(all.equal(
  switch(as.character(family),
         lognormal   = sigma_dgm^2,
         loglogistic = sigma_dgm^2 * pi^2 / 3,
         cox         = pi^2 / (6 * rho_dgm^2)),
  sigma^2)))

# tau is the shared AFT-scale frailty SD target.  u_AFT = -u_Cox/rho, so the Cox
# arm is simulated at rho_dgm * tau to give an implied AFT-scale SD of exactly tau.
tau_dgm <- if (as.character(family) == "cox") rho_dgm * tau else tau

# AFT-scale targets under the Cox DGM.  Inverting the Weibull baseline gives
#   log T = -(log lambda_0)/rho - eta_Cox/rho - u/rho + (log E)/rho,  E ~ Exp(1)
# so every effect maps to log time by division by -rho.  The intercept also
# absorbs the baseline scale -(log lambda_0)/rho and, because the fitted
# Gaussian error is mean zero while (log E)/rho has mean -gamma_EM/rho, that
# offset too.
EULER_MASCHERONI <- -digamma(1)          # 0.5772157

aft_beta_target <- function(beta_cox, rho) -beta_cox / rho

aft_gamma_target <- function(gamma_cox, rho, lambda_0) {
  c(-(log(lambda_0) + gamma_cox[1] + EULER_MASCHERONI) / rho,
    -gamma_cox[-1] / rho)
}

###############################################################
## run simulations
###############################################################

# collect one-row data frames per successful iteration
coef_list <- vector("list", length = N_iter)
info_list <- vector("list", length = N_iter)

###############################################################
## Marginal test set: one draw shared by every iteration.
## New clusters -> frailty = 0 at prediction time.
###############################################################
set.seed(scenario)

# Deliberately NOT a mirror of the training design.  Predictions here set the
# frailty to 0, so the test clusters exist only to inject between-cluster
# variance into the test outcomes; since this set is drawn ONCE and reused by
# every iteration, few clusters would mean a small sample of u_j whose
# idiosyncrasy biases the whole scenario and never averages out.  Many small
# clusters (~200 of ~5) fix that at no cost.
n_cluster_test <- max(50L, N_test %/% 5L)
if (family %in% c("lognormal", "loglogistic")) {
  sim_test <- simulate_AFT_interaction(
    family       = as.character(family),
    N_total      = N_test,
    n_cluster    = n_cluster_test,
    nS           = nS,
    beta_type    = as.character(beta_type),
    gamma        = gamma_dgm,
    tau          = tau,
    sigma        = sigma_dgm,
    censor_rate  = censor_rate
  )
} else {
  sim_test <- simulate_Cox_interaction(
    N_total      = N_test,
    n_cluster    = n_cluster_test,
    nS           = nS,
    beta_type    = as.character(beta_type),
    gamma        = gamma_dgm,
    lambda_0     = 0.1,
    rho          = rho_dgm,
    beta_scale   = rho_dgm,   # implied AFT surface matches the AFT arms
    tau          = tau_dgm,   # = rho_dgm*tau: implied AFT-scale SD is tau
    censor_rate  = censor_rate
  )
}
test_data <- sim_test$data
Z_test    <- model.matrix(~ Z1 + Z2, data = test_data)
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

# Precompute survival objects used by UnoC (fixed across all iterations)
Surv_test <- survival::Surv(test_data$Y, test_data$delta)

# ---------------------------------------------------------------------------
# Precompute the surface evaluation grid once (constant within a scenario).
# s_grid and age_range depend only on scenario parameters, not random seed.
# beta_true_grid is also constant (same true surface for all iterations).
# ---------------------------------------------------------------------------
n_eval_s   <- 30
n_eval_age <- 30
age_range  <- sim_test$age_range
s_eval     <- seq(0, 1, length.out = n_eval_s)          # tmax = 1 always
age_eval   <- seq(age_range[1], age_range[2], length.out = n_eval_age)

grid_df        <- expand.grid(s_grid = s_eval, age = age_eval)
a_grid_scaled  <- (grid_df$age - age_range[1]) / diff(age_range)
beta_true_dgm  <- mapply(sim_test$beta_fn, grid_df$s_grid, a_grid_scaled)

# AFT-scale target surface: unchanged for the AFT families, -beta_Cox/rho for
# the Cox DGM, whose beta_fn lives on the log-hazard scale.  Without this the
# Cox surface metrics measure the scale gap rather than estimation error
# (mean(beta_additive) = 1 and rho = 1.5 give an offset of 1 + 1/1.5 = 1.67,
# which is exactly the bias previously reported for the Cox arm).
beta_true_grid <- if (as.character(family) == "cox") {
  aft_beta_target(beta_true_dgm, rho_dgm)
} else {
  beta_true_dgm
}

n_point  <- nrow(grid_df)
mat_grid <- function() matrix(NA_real_, N_iter, n_point)
mat_s    <- function() matrix(NA_real_, N_iter, n_eval_s)
surf     <- list(mean = mat_grid(), sd = mat_grid(),
                 q025 = mat_grid(), q975 = mat_grid())
# ni_surf  <- list(mean = mat_s(), sd = mat_s(), q025 = mat_s(), q975 = mat_s())

# Per-iteration seeds
seeds <- sample.int(1e8, N_iter)

for(iter in 1:N_iter){
  cat("Scenario", scenario, "| Iteration", iter, "\n")
  # set seed
  set.seed(seeds[[iter]])
  
  # simulate training data — retry up to 10 times if cluster-size
  # sampling fails (should be prevented by the repeat loop in the
  # simulate functions, but kept as a belt-and-suspenders guard)
  sim_data <- NULL
  for (.attempt in seq_len(10L)) {
    sim_data <- tryCatch({
      if (family %in% c("lognormal", "loglogistic")) {
        simulate_AFT_interaction(
          family       = as.character(family),
          N_total      = N_total,
          n_cluster    = n_cluster,
          nS           = nS,
          beta_type    = as.character(beta_type),
          gamma        = gamma_dgm,
          tau          = tau,
          sigma        = sigma_dgm,
          censor_rate  = censor_rate
        )
      } else {
        simulate_Cox_interaction(
          N_total      = N_total,
          n_cluster    = n_cluster,
          nS           = nS,
          beta_type    = as.character(beta_type),
          gamma        = gamma_dgm,
          lambda_0     = 0.1,
          rho          = rho_dgm,
          beta_scale   = rho_dgm,   # implied AFT surface matches the AFT arms
          tau          = tau_dgm,   # = rho_dgm*tau: implied AFT-scale SD is tau
          censor_rate  = censor_rate
        )
      }
    }, error = function(e) {
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
    ## fit functional frailty AFT  model
    ###############################################################
    
    # extract elements from simulated data
    data <- sim_data$data
    s_grid = sim_data$s_grid
    Z <- model.matrix(~ Z1 + Z2, data = data)
    
    # run Gibbs sampler
    tic()
    fit <- gibbs_frailty_interaction(
      time = data$Y,
      status = data$delta,
      cluster_id = data$cluster_id,
      Z = Z,
      X = data$X,
      age = data$age,
      s_grid = s_grid,
      # tuning / priors
      K_s = K_s_fit,
      K_age = K_age_fit,
      bs_s = bs_s_fit,
      lambda_init = 1,         # smoothing parameter
      A_lambda = 1, B_lambda = 1,
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
      A_sigma2 = 3, B_sigma2 = 2,  # IG(A,B) for sigma^2 (shape A, scale B)
      # MCMC
      n_iter = n_iter_fit,
      n_burn = n_burn_fit,
      n_thin = 1,
      verbose = TRUE
    )
    
    time_stamp <- toc(quiet = TRUE)
    time <- time_stamp$toc - time_stamp$tic
    
    ###############################################################
    ## reconstruct posterior surface beta(s, a) on the precomputed grid
    ## fit$b holds S x K_total draws; the surface is b %*% t(Phi_grid)
    ## where Phi_grid = PredictMat(sm_tensor, grid_df)
    ## grid_df, a_grid_scaled, beta_true_grid are constant across iterations
    ## (precomputed outside the loop); only Phi_grid changes with each fit.
    ###############################################################
    
    # Basis matrix on the eval grid (reusing the sm_tensor stored in the fit)
    Phi_grid   <- PredictMat(fit$meta$sm_tensor, grid_df)  # n_grid x K_total
    
    # Posterior draws of the surface: S x n_grid
    surface_draws <- fit$b %*% t(Phi_grid)
    
    surface_mean  <- colMeans(surface_draws)
    surface_sd    <- apply(surface_draws, 2, sd)
    surface_q025  <- apply(surface_draws, 2, quantile, 0.025, na.rm = TRUE)
    surface_q975  <- apply(surface_draws, 2, quantile, 0.975, na.rm = TRUE)

    # Simultaneous (max-statistic) band: beta_hat +/- qd * sd(t).  Pointwise
    # coverage is not a statement about the surface -- the pointwise band holds
    # the ENTIRE true surface only ~58% of the time despite ~0.99 per point.
    cma            <- cma_band(surface_draws, level = 0.95)
    beta_cover_sim <- all(beta_true_grid >= cma$lower & beta_true_grid <= cma$upper)
    beta_cma_qd    <- cma$qd
    
    ###############################################################
    ## surface estimation metrics
    ###############################################################
    
    # The three surfaces carry different signal: sd(int X beta) is 1.06 for
    # "additive" but 0.61 for "bilinear" and "ridge", so a raw ISE comparison
    # across beta_type confounds "harder surface" with "weaker signal".
    # beta_ise_rel divides by mean(beta_true^2) and is scale-free.
    beta_true_ss   <- mean(beta_true_grid^2)
    beta_resid     <- surface_mean - beta_true_grid
    beta_ise       <- mean(beta_resid^2)           # approx integral over the 2D grid
    beta_ise_rel   <- beta_ise / beta_true_ss
    beta_bias_mean <- mean(beta_resid)             # signed mean bias
    beta_bias_max  <- max(abs(beta_resid))         # max absolute pointwise bias
    beta_cover_pw  <- mean(                        # pointwise 95% CrI coverage
      (beta_true_grid >= surface_q025) & (beta_true_grid <= surface_q975)
    )
    
    ###############################################################
    ## scalar parameter metrics
    ###############################################################
    
    # --- gamma (intercept + Z1 + Z2) ---
    # gamma_true = DGM scale; gamma_implied = AFT target, and what bias/SE/
    # coverage are computed against.  They differ only for the Cox DGM.
    gamma_true    <- sim_data$gamma
    gamma_implied <- if (as.character(family) == "cox") {
      aft_gamma_target(gamma_true, sim_data$rho, sim_data$lambda_0)
    } else {
      gamma_true    # AFT DGM: same scale, no correction needed
    }
    gamma_est   <- colMeans(fit$gamma)
    gamma_bias  <- gamma_est  - gamma_implied
    gamma_se    <- (gamma_est - gamma_implied)^2
    q_gamma     <- apply(fit$gamma, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    gamma_cover <- (gamma_implied >= q_gamma[1, ]) & (gamma_implied <= q_gamma[2, ])
    
    # --- tau2 ---
    # tau2_implied = tau2_true/rho^2 for Cox, which given tau_dgm = rho*tau
    # evaluates to tau^2 -- the same target as the AFT arms.
    tau2_true    <- sim_data$tau^2
    tau2_implied <- if (as.character(family) == "cox") {
      tau2_true / sim_data$rho^2
    } else {
      tau2_true     # AFT DGM: same scale, no correction needed
    }
    # Posterior MEDIAN, not mean: the tau2 posterior is right-skewed (skewness
    # 1.2-1.4 at J = 25), so the mean sits above the truth by construction.  At
    # tau = 1, J = 25 the median gives 1.119 against a mean of 1.188.
    tau2_est   <- median(fit$tau2)
    tau2_bias  <- tau2_est - tau2_implied   # compare against implied target
    tau2_se    <- (tau2_est - tau2_implied)^2
    tau2_cover <- (tau2_implied >= quantile(fit$tau2, 0.025)) &
      (tau2_implied <= quantile(fit$tau2, 0.975))
    
    # --- sigma2 ---
    # The target is the log-time RESIDUAL variance, which is what the Gaussian
    # fit estimates (its full conditional targets the residual second moment),
    # not the DGM's scale parameter.  All three evaluate to sigma^2 by the
    # sigma_dgm / rho_dgm construction, asserted outside the loop.
    sigma2_true  <- switch(
      as.character(family),
      lognormal   = sim_data$sigma^2,
      loglogistic = sim_data$sigma^2 * pi^2 / 3,
      cox         = pi^2 / (6 * sim_data$rho^2)
    )
    sigma2_est   <- mean(fit$sigma2)
    sigma2_bias  <- sigma2_est - sigma2_true
    sigma2_se    <- (sigma2_est - sigma2_true)^2
    sigma2_cover <- (sigma2_true >= quantile(fit$sigma2, 0.025)) &
      (sigma2_true <= quantile(fit$sigma2, 0.975))
    
    ###############################################################
    ## out-of-sample predictive metrics on test set
    ## New clusters -> cluster_id_new = NULL (marginal, frailty = 0)
    ###############################################################
    
    # --- survival probability draws on test set ---
    pred_surv <- predict_gibbs_interaction(
      fit            = fit,
      X_new          = test_data$X,
      age_new        = test_data$age,
      Z_new          = Z_test,
      cluster_id_new = NULL,
      type           = "survival",
      times          = tgrid,
      level          = 0.95,
      quantiles      = FALSE      # only $mean is used; skip the credible bands
    )
    
    # --- risk score for the C-index ---
    # Negative posterior-mean linear predictor.  Rank-equivalent to
    # -rowSums(S) for an AFT model but cannot underflow: a saturated survival
    # matrix collapses -rowSums() to one repeated value, and UnoC scores only
    # strict lp_i > lp_j, so ties return C = 0 rather than 0.5.
    pred_link  <- predict_gibbs_interaction(
      fit            = fit,
      X_new          = test_data$X,
      age_new        = test_data$age,
      Z_new          = Z_test,
      cluster_id_new = NULL,
      type           = "link",
      level          = 0.95,
      quantiles      = FALSE
    )
    risk_score <- -pred_link$mean            # length N_test; higher = higher risk
    
    # Uno's IPCW C-index.
    # Upweights comparable pairs by the inverse probability of censoring,
    # making it consistent even under dependent censoring.
    # UnoC(Surv.rsp, Surv.rsp.new, lpnew) where lpnew = risk score for test set.
    Surv_train <- survival::Surv(data$Y, data$delta)  # training set (changes each iter)
    c_index    <- survAUC::UnoC(
      Surv.rsp     = Surv_train,
      Surv.rsp.new = Surv_test,
      lpnew        = risk_score,
      time         = tau_eval   # explicit truncation; default would be max(test time)
    )
    
    # IPCW Integrated Brier Score.
    # tgrid spans 5th-95th percentile of test EVENT times (precomputed outside
    # the loop) to include tails where lognormal and Weibull diverge most.
    ibs <- as.numeric(cal_IPCW_Brier(
      S_mat       = pred_surv$mean,
      time_test   = test_data$Y,
      event_test  = test_data$delta,
      time_train  = data$Y,
      event_train = data$delta,
      tgrid       = tgrid,
      G_data      = "test",   # censoring KM from the sample being weighted
      eps_G       = eps_G
    ))
    
    ###############################################################
    ## out-of-sample predictive metrics — CONDITIONAL (known clusters)
    ##
    ## New subjects in the SAME clusters as the training data, sharing its u_j
    ## via the DGM's `u` argument and held out of the fit.  Scored twice with the
    ## same fit and subjects: *_known uses u_hat_j, *_known_marg sets frailty to
    ## 0; the difference is the frailty's predictive contribution.  No leakage —
    ## these subjects' outcomes never entered the sampler.  Unlike the marginal
    ## test set this one cannot be fixed across iterations, so expect more noise.
    ###############################################################

    sim_known <- if (family %in% c("lognormal", "loglogistic")) {
      simulate_AFT_interaction(
        family = as.character(family), N_total = N_known, n_cluster = n_cluster,
        nS = nS, beta_type = as.character(beta_type), gamma = gamma_dgm,
        tau = tau, sigma = sigma_dgm, censor_rate = censor_rate,
        u = sim_data$u                 # share the training frailties
      )
    } else {
      simulate_Cox_interaction(
        N_total = N_known, n_cluster = n_cluster, nS = nS,
        beta_type = as.character(beta_type), gamma = gamma_dgm,
        lambda_0 = 0.1, rho = rho_dgm, tau = tau_dgm, censor_rate = censor_rate,
        beta_scale = rho_dgm,   # implied AFT surface = -beta_fn, matching the AFT arms
        u = sim_data$u                 # share the training frailties
      )
    }
    known_data <- sim_known$data
    Z_known    <- model.matrix(~ Z1 + Z2, data = known_data)
    Surv_known <- survival::Surv(known_data$Y, known_data$delta)

    # predict_gibbs_interaction() silently falls back to frailty = 0 for any
    # cluster id it cannot match, which would relabel marginal predictions as
    # conditional.  Fail loudly instead.
    stopifnot(all(as.character(known_data$cluster_id) %in% fit$meta$cluster_levels))

    # Reuse the marginal horizon so the two evaluations are on a common scale,
    # trimmed to where the known-cluster set's own censoring KM supports IPCW.
    km_known    <- survival::survfit(survival::Surv(known_data$Y, 1 - known_data$delta) ~ 1)
    G_known     <- stats::stepfun(km_known$time, c(1, km_known$surv))
    tgrid_known <- tgrid[G_known(tgrid) >= eps_G]
    stopifnot(length(tgrid_known) >= 2L)
    tau_known   <- max(tgrid_known)

    eval_known <- function(cid) {
      lk <- predict_gibbs_interaction(
        fit = fit, X_new = known_data$X, age_new = known_data$age,
        Z_new = Z_known, cluster_id_new = cid, type = "link", level = 0.95,
        quantiles = FALSE)
      sv <- predict_gibbs_interaction(
        fit = fit, X_new = known_data$X, age_new = known_data$age,
        Z_new = Z_known, cluster_id_new = cid, type = "survival",
        times = tgrid_known, level = 0.95, quantiles = FALSE)
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
    # the orders agree.  For the Cox DGM u is on the log-hazard scale and
    # u_AFT = -u_Cox/rho (see the sign-flip note in the audit), so put the truth
    # on the AFT scale before correlating.
    u_hat  <- colMeans(fit$u)
    u_true <- sim_data$u[as.integer(fit$meta$cluster_levels)]
    if (as.character(family) == "cox") u_true <- -u_true / sim_data$rho
    cor_u  <- suppressWarnings(cor(u_hat, u_true))

    # ###############################################################
    # ## COMPARATOR — no-interaction functional frailty AFT
    # ##
    # ## Identical to the fitted model except the coefficient is beta(s) rather
    # ## than beta(s, age): same error, frailty, augmentation, priors, MCMC length
    # ## and s-basis, fitted to the same data and scored on the same grid and test
    # ## sets, so every ni_* column pairs with its unprefixed counterpart.
    # ## beta_hat(s) is broadcast along age -- the restriction the comparator
    # ## imposes -- so ni_beta_ise vs beta_ise is the cost of ignoring the
    # ## interaction.
    # ###############################################################
#
    # tic()
    # fit_ni <- gibbs_functional_frailty(
      # time = data$Y, status = data$delta, cluster_id = data$cluster_id,
      # Z = Z, X = data$X, s_grid = s_grid,
      # K = K_s_fit, basis_type = bs_s_fit,
      # lambda_init = 1, A_lambda = 1, B_lambda = 1,
      # var_gamma = 100,
      # A_tau2 = 0.01, B_tau2 = 0.01,
      # A_sigma2 = 3, B_sigma2 = 2,
      # n_iter = n_iter_fit, n_burn = n_burn_fit, n_thin = 1, verbose = FALSE
    # )
    # ts_ni   <- toc(quiet = TRUE)
    # time_ni <- ts_ni$toc - ts_ni$tic
#
    # # --- beta_hat(s) on the evaluation grid, held constant across age ---
    # # Rebuild the sampler's own smooth on s_grid, then evaluate at s_eval.
    # sm_ni <- smoothCon(s(s_grid, k = K_s_fit, bs = bs_s_fit),
                       # data = data.frame(s_grid = s_grid), absorb.cons = FALSE)[[1]]
    # Phi_ni <- PredictMat(sm_ni, data.frame(s_grid = s_eval))   # n_eval_s x K
    # stopifnot(ncol(Phi_ni) == ncol(fit_ni$b))
#
    # beta_ni_draws <- fit_ni$b %*% t(Phi_ni)                    # S x n_eval_s
    # ni_beta_mean  <- colMeans(beta_ni_draws)
    # ni_beta_q025  <- apply(beta_ni_draws, 2, quantile, 0.025, na.rm = TRUE)
    # ni_beta_q975  <- apply(beta_ni_draws, 2, quantile, 0.975, na.rm = TRUE)
#
    # # Band on the n_eval_s distinct columns, not the broadcast grid: duplicated
    # # columns change neither a per-column mean/sd nor a maximum, so qd and the
    # # band are identical without allocating the S x n_grid matrix.
    # cma_ni            <- cma_band(beta_ni_draws, level = 0.95)
    # ni_beta_sd        <- apply(beta_ni_draws, 2, sd)
#
    # # Broadcast along age (summaries only — never the S x n_grid draw matrix)
    # s_idx           <- match(grid_df$s_grid, s_eval)
    # ni_surface_mean <- ni_beta_mean[s_idx]
    # ni_surface_q025 <- ni_beta_q025[s_idx]
    # ni_surface_q975 <- ni_beta_q975[s_idx]
#
    # ni_beta_cover_sim <- all(beta_true_grid >= cma_ni$lower[s_idx] &
                             # beta_true_grid <= cma_ni$upper[s_idx])
    # ni_beta_cma_qd    <- cma_ni$qd
#
    # ni_resid          <- ni_surface_mean - beta_true_grid
    # ni_beta_ise       <- mean(ni_resid^2)
    # ni_beta_ise_rel   <- ni_beta_ise / beta_true_ss
    # ni_beta_bias_mean <- mean(ni_resid)
    # ni_beta_bias_max  <- max(abs(ni_resid))
    # ni_beta_cover_pw  <- mean((beta_true_grid >= ni_surface_q025) &
                              # (beta_true_grid <= ni_surface_q975))
#
    # # --- scalar parameters, against the same targets ---
    # ni_gamma_est   <- colMeans(fit_ni$gamma)
    # ni_gamma_bias  <- ni_gamma_est - gamma_implied
    # ni_gamma_se    <- (ni_gamma_est - gamma_implied)^2
    # q_ni_gamma     <- apply(fit_ni$gamma, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    # ni_gamma_cover <- (gamma_implied >= q_ni_gamma[1, ]) & (gamma_implied <= q_ni_gamma[2, ])
#
    # ni_tau2_est   <- median(fit_ni$tau2)
    # ni_tau2_bias  <- ni_tau2_est - tau2_implied
    # ni_tau2_se    <- (ni_tau2_est - tau2_implied)^2
    # ni_tau2_cover <- (tau2_implied >= quantile(fit_ni$tau2, 0.025)) &
      # (tau2_implied <= quantile(fit_ni$tau2, 0.975))
#
    # ni_sigma2_est   <- mean(fit_ni$sigma2)
    # ni_sigma2_bias  <- ni_sigma2_est - sigma2_true
    # ni_sigma2_se    <- (ni_sigma2_est - sigma2_true)^2
    # ni_sigma2_cover <- (sigma2_true >= quantile(fit_ni$sigma2, 0.025)) &
      # (sigma2_true <= quantile(fit_ni$sigma2, 0.975))
#
    # # --- predictive: same test sets, same horizons, same risk-score convention ---
    # ni_pred_marg <- {
      # lk <- predict_gibbs_frailty(fit_ni, X_new = test_data$X, Z_new = Z_test,
              # cluster_id_new = NULL, type = "link", level = 0.95, quantiles = FALSE)
      # sv <- predict_gibbs_frailty(fit_ni, X_new = test_data$X, Z_new = Z_test,
              # cluster_id_new = NULL, type = "survival", times = tgrid,
              # level = 0.95, quantiles = FALSE)
      # c(c_index = survAUC::UnoC(Surv.rsp = Surv_train, Surv.rsp.new = Surv_test,
                                # lpnew = -lk$mean, time = tau_eval),
        # ibs = as.numeric(cal_IPCW_Brier(sv$mean, test_data$Y, test_data$delta,
                                        # data$Y, data$delta, tgrid,
                                        # G_data = "test", eps_G = eps_G)))
    # }
#
    # eval_known_ni <- function(cid) {
      # lk <- predict_gibbs_frailty(fit_ni, X_new = known_data$X, Z_new = Z_known,
              # cluster_id_new = cid, type = "link", level = 0.95, quantiles = FALSE)
      # sv <- predict_gibbs_frailty(fit_ni, X_new = known_data$X, Z_new = Z_known,
              # cluster_id_new = cid, type = "survival", times = tgrid_known,
              # level = 0.95, quantiles = FALSE)
      # c(c_index = survAUC::UnoC(Surv.rsp = Surv_train, Surv.rsp.new = Surv_known,
                                # lpnew = -lk$mean, time = tau_known),
        # ibs = as.numeric(cal_IPCW_Brier(sv$mean, known_data$Y, known_data$delta,
                                        # data$Y, data$delta, tgrid_known,
                                        # G_data = "test", eps_G = eps_G)))
    # }
    # ni_cond <- eval_known_ni(known_data$cluster_id)
    # ni_marg <- eval_known_ni(NULL)

    ###############################################################
    ## assemble surface data frame  (n_grid rows, one per eval point)
    ###############################################################
    
    # # Rows for the surface matrices.  The comparator is stored at its native
    # # n_eval_s resolution, not broadcast along age.
    surface_row <- list(mean = surface_mean, sd = surface_sd,
                        q025 = surface_q025, q975 = surface_q975)
    # ni_row      <- list(mean = ni_beta_mean, sd = ni_beta_sd,
                        # q025 = ni_beta_q025, q975 = ni_beta_q975)
    
    ###############################################################
    ## assemble one-row result data frames
    ###############################################################
    
    # ---- coefficient / estimation summary ----
    df_coef <- data.frame(
      # surface, scored against the AFT-scale target
      beta_ise       = beta_ise,
      beta_ise_rel   = beta_ise_rel,   # ISE / mean(beta_true^2); comparable across beta_type
      beta_bias_mean = beta_bias_mean,
      beta_bias_max  = beta_bias_max,
      beta_cover_pw  = beta_cover_pw,
      # simultaneous: is the WHOLE true surface inside the band? (nominal 0.95)
      beta_cover_sim = beta_cover_sim,
      beta_cma_qd    = beta_cma_qd,    # critical value; pointwise would use 1.96
      # frailty variance; tau2_implied is the AFT-scale target for bias/coverage
      tau2_true    = tau2_true,
      tau2_implied = tau2_implied,
      tau2_est     = tau2_est,
      tau2_bias    = tau2_bias,       # tau2_est - tau2_implied
      tau2_se      = tau2_se,
      tau2_cover   = tau2_cover,      # I(tau2_implied in 95% CrI of fit$tau2)
      # residual variance; = sigma^2 for every family by construction
      sigma2_true  = sigma2_true,
      sigma2_est   = sigma2_est,
      sigma2_bias  = sigma2_bias,
      sigma2_se    = sigma2_se,
      sigma2_cover = sigma2_cover,
      # predictive: discrimination (C-index) and calibration (IBS)
      c_index      = c_index,
      ibs          = ibs,
      # predictive — conditional, known clusters (see the block above).
      # *_known uses u_hat_j; *_known_marg is the same subjects with frailty 0.
      # The difference is the frailty's predictive contribution.
      c_index_known      = m_cond[["c_index"]],
      ibs_known          = m_cond[["ibs"]],
      c_index_known_marg = m_marg[["c_index"]],
      ibs_known_marg     = m_marg[["ibs"]],
      # frailty recovery: cor(u_hat_j, u_j) on the AFT scale
      cor_u              = cor_u,
      # diagnostics: a collapsed risk score drives Uno's C to 0 via ties, so
      # record enough to spot it instead of silently averaging it in
      risk_sd       = sd(risk_score),
      risk_n_unique = length(unique(risk_score))
      # # ---- no-interaction comparator, beta(s) instead of beta(s, age) ----
      # # Same data, same priors, same grid, same test sets and same targets, so
      # # each ni_* pairs directly with the column of the same name above.
      # ni_beta_ise       = ni_beta_ise,
      # ni_beta_ise_rel   = ni_beta_ise_rel,
      # ni_beta_bias_mean = ni_beta_bias_mean,
      # ni_beta_bias_max  = ni_beta_bias_max,
      # ni_beta_cover_pw  = ni_beta_cover_pw,
      # ni_beta_cover_sim = ni_beta_cover_sim,
      # ni_beta_cma_qd    = ni_beta_cma_qd,
      # ni_tau2_est       = ni_tau2_est,
      # ni_tau2_bias      = ni_tau2_bias,
      # ni_tau2_se        = ni_tau2_se,
      # ni_tau2_cover     = ni_tau2_cover,
      # ni_sigma2_est     = ni_sigma2_est,
      # ni_sigma2_bias    = ni_sigma2_bias,
      # ni_sigma2_se      = ni_sigma2_se,
      # ni_sigma2_cover   = ni_sigma2_cover,
      # ni_c_index             = ni_pred_marg[["c_index"]],
      # ni_ibs                 = ni_pred_marg[["ibs"]],
      # ni_c_index_known       = ni_cond[["c_index"]],
      # ni_ibs_known           = ni_cond[["ibs"]],
      # ni_c_index_known_marg  = ni_marg[["c_index"]],
      # ni_ibs_known_marg      = ni_marg[["ibs"]]
    )
    
    # gamma columns added dynamically (intercept = index 0, Z1 = 1, Z2 = 2)
    # gamma_true_k    : DGM scale (log-hazard for Cox, log-time for AFT)
    # gamma_implied_k : AFT-scale target; bias/se/cover are relative to this
    idx <- seq_along(gamma_est) - 1L
    df_coef[paste0("gamma_true_",    idx)] <- as.list(gamma_true)
    df_coef[paste0("gamma_implied_", idx)] <- as.list(gamma_implied)
    df_coef[paste0("gamma_est_",     idx)] <- as.list(gamma_est)
    df_coef[paste0("gamma_bias_",    idx)] <- as.list(gamma_bias)
    df_coef[paste0("gamma_se_",      idx)] <- as.list(gamma_se)
    df_coef[paste0("gamma_cover_",   idx)] <- as.list(gamma_cover)
    # df_coef[paste0("ni_gamma_est_",   idx)] <- as.list(ni_gamma_est)
    # df_coef[paste0("ni_gamma_bias_",  idx)] <- as.list(ni_gamma_bias)
    # df_coef[paste0("ni_gamma_se_",    idx)] <- as.list(ni_gamma_se)
    # df_coef[paste0("ni_gamma_cover_", idx)] <- as.list(ni_gamma_cover)
    
    # ---- run-level info ----
    df_info <- data.frame(
      scenario    = scenario,
      iter        = iter,
      seed        = seeds[[iter]],
      family      = as.character(family),
      N_total     = N_total,
      n_cluster   = n_cluster,
      nS          = nS,
      beta_type   = as.character(beta_type),
      tau         = tau,      # AFT-scale target, common across families
      tau_dgm     = tau_dgm,  # SD passed to the DGM (= rho * tau for Cox)
      sigma       = sigma,      # log-time residual SD target, common to all
      sigma_dgm   = sigma_dgm,  # scale param handed to the AFT DGM (NA for Cox)
      rho         = rho_dgm,    # Weibull shape used by the Cox DGM (NA otherwise)
      censor_rate = censor_rate,
      event_rate  = mean(data$delta),
      time        = time # Gibbs seconds, interaction model
      # time_ni     = time_ni    # Gibbs seconds, no-interaction comparator
    )
    
    list(info = df_info, coef = df_coef,
         surface_row = surface_row)
    
  }, error = function(e) {
    warning(sprintf(
      "Iteration %d skipped due to error:\n  %s",
      iter, e$message
    ))
    NULL
  })
  
  if (!is.null(res)) {
    info_list[[iter]] <- res$info
    coef_list[[iter]] <- res$coef
    for (k in names(surf))    surf[[k]][iter, ]    <- res$surface_row[[k]]
    # for (k in names(ni_surf)) ni_surf[[k]][iter, ] <- res$ni_row[[k]]
  }
} # end for loop

# drop NULL entries and bind into data frames
# Row i of every matrix is iteration i; use result$coef$iter to select the
# iterations that completed.  grid describes the columns of surface$*, and
# s_eval the columns of ni_surface$* (whose surface is constant in age).
result <- list(
  info       = dplyr::bind_rows(Filter(Negate(is.null), info_list)),
  coef       = dplyr::bind_rows(Filter(Negate(is.null), coef_list)),
  grid       = data.frame(s          = grid_df$s_grid,
                          age        = grid_df$age,
                          age_scaled = a_grid_scaled,
                          beta_true  = beta_true_grid),
  s_eval     = s_eval,
  surface    = surf
  # ni_surface = ni_surf
)

###############################################################
## save result
###############################################################
# Suffix the folder with the study name.  Both runners are date-stamped, so
# without this they would write outputs/<date>/<scenario>.RDA to the SAME path
# and whichever finished second would silently overwrite the other.
Date   <- gsub("-", "", Sys.Date())
outdir <- file.path(here::here("outputs"), paste0(Date, "_interaction"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

filename <- file.path(outdir, paste0(scenario, ".RDA"))
save(result, file = filename)

###############################################################
## end sim
###############################################################