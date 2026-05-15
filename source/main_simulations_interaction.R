####################################################################
# Weijia Qian
# Jan 10, 2026
#
# Simulate survival data under a functional frailty AFT or Cox model
# with a bivariate surface beta(s, age), fit gibbs_frailty_interaction(),
# and extract estimation + predictive metrics.
#
# Metrics saved per iteration
# ---------------------------
# Estimation — surface beta(s, a):
#   beta_ise         integrated squared error over the eval grid
#   beta_bias_mean   mean pointwise bias
#   beta_bias_max    max absolute pointwise bias
#   beta_cover_pw    pointwise 95% CrI coverage proportion
#
# Estimation — scalar parameters:
#   gamma_{est,bias,se,cover}_k   for each gamma_k (k = 0,1,2)
#     NOTE: for Cox DGM, gamma is on the log-hazard scale in the DGM but
#     log-time scale in the fitted AFT model.  Comparison is approximate.
#   tau2_{true,implied,est,bias,se,cover}
#     For Cox DGM, tau2_implied = tau2_true / rho^2 (AFT-equivalent scale).
#     Bias/coverage are computed against tau2_implied.
#   sigma2_{true,est,bias,se,cover}
#     For Cox DGM, sigma2_true = NA (no direct analogue).
#
# Predictive — out-of-sample test set (N_test = 1000 new subjects,
#              new clusters, same DGM parameters, marginal frailty = 0):
#   c_index          Uno's IPCW C-index (survAUC::UnoC)
#   ibs              IPCW integrated Brier score (cal_Brier.R / cal_IPCW_Brier)
#
# Computational:
#   time             wall-clock seconds for the Gibbs sampler
####################################################################

suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(refund))
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
source(here("source", "cal_Brier.R"))
source(here("source", "simulate_AFT_interaction.R"))
source(here("source", "simulate_Cox_interaction.R"))

###############################################################
## set simulation design elements
###############################################################
family = c("lognormal", "weibull", "cox")
N_total = c(4000)
n_cluster = c(25, 50, 75)
nS = c(100)
beta_type = c("additive")
# beta_type = c("additive", "bilinear", "ridge")
tau = c(0.2, 0.5, 1)
sigma = c(1)       # AFT error SD (= 1/rho for Weibull; ignored for Cox DGM)
rho   = c(1.5)     # Weibull shape for Cox DGM (ignored for AFT DGMs)
censor_rate = c(0.5, 0.75)
N_iter = 200
N_test = 1000      # out-of-sample test set size (new subjects, new clusters)

params = expand.grid(family = family,
                     N_total = N_total,
                     n_cluster = n_cluster,
                     nS = nS,
                     beta_type = beta_type,
                     tau = tau,
                     sigma = sigma,
                     rho = rho,
                     censor_rate = censor_rate)

## define number of simulations and parameter scenarios
if(doLocal) {
  scenario = 30
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
rho         = params$rho[scenario]
censor_rate = params$censor_rate[scenario]

###############################################################
## run simulations
###############################################################

# collect one-row data frames per successful iteration
coef_list    <- vector("list", length = N_iter)
info_list    <- vector("list", length = N_iter)
surface_list <- vector("list", length = N_iter)  # n_grid-row df per iter
seeds <- sample.int(1e8, N_iter)

###############################################################
## generate a single test set (shared across all iterations)
## Fixed seed derived from scenario so results are reproducible
## across runs. New clusters -> frailty = 0 at prediction time.
###############################################################
set.seed(scenario * 1000L + 99L)
n_cluster_test <- max(5L, round(N_test / (N_total / n_cluster)))
if (family %in% c("lognormal", "weibull")) {
  sim_test <- simulate_AFT_interaction(
    family       = as.character(family),
    N_total      = N_test,
    n_cluster    = n_cluster_test,
    nS           = nS,
    beta_type    = as.character(beta_type),
    gamma        = c(0.5, 0.3, -0.2),
    tau          = tau,
    sigma        = sigma,
    censor_rate  = censor_rate
  )
} else {
  sim_test <- simulate_Cox_interaction(
    N_total      = N_test,
    n_cluster    = n_cluster_test,
    nS           = nS,
    beta_type    = as.character(beta_type),
    gamma        = c(0.5, 0.3, -0.2),
    lambda_0     = 0.1,
    rho          = rho,
    tau          = tau,
    censor_rate  = censor_rate
  )
}
test_data <- sim_test$data
Z_test    <- model.matrix(~ Z1 + Z2, data = test_data)
cat("Test set generated: N =", nrow(test_data),
    "| event rate =", round(mean(test_data$delta), 3), "\n")

# ---------------------------------------------------------------------------
# Precompute tgrid once (test_data is fixed across all iterations).
# Use quantiles of EVENT times only (not censored) so the grid spans the
# actual event time distribution, including the tails where lognormal and
# Weibull differ most.  5th-95th percentile of event times to capture tails.
# ---------------------------------------------------------------------------
event_times_test <- test_data$Y[test_data$delta == 1]
tgrid <- quantile(event_times_test, probs = seq(0.05, 0.95, by = 0.05))
tgrid <- sort(unique(tgrid))

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
a_eval     <- (age_eval - age_range[1]) / diff(age_range)

grid_df        <- expand.grid(s_grid = s_eval, age = age_eval)
a_grid_scaled  <- (grid_df$age - age_range[1]) / diff(age_range)
beta_true_grid <- mapply(sim_test$beta_fn, grid_df$s_grid, a_grid_scaled)

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
      if (family %in% c("lognormal", "weibull")) {
        simulate_AFT_interaction(
          family       = as.character(family),
          N_total      = N_total,
          n_cluster    = n_cluster,
          nS           = nS,
          beta_type    = as.character(beta_type),
          gamma        = c(0.5, 0.3, -0.2),
          tau          = tau,
          sigma        = sigma,
          censor_rate  = censor_rate
        )
      } else {
        simulate_Cox_interaction(
          N_total      = N_total,
          n_cluster    = n_cluster,
          nS           = nS,
          beta_type    = as.character(beta_type),
          gamma        = c(0.5, 0.3, -0.2),
          lambda_0     = 0.1,
          rho          = rho,
          tau          = tau,
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
      K_s = 15,
      K_age = 10,
      bs_s = "cr",    # "cc" for real NHANES data, "cr" for simulation
      lambda_init = 1,         # smoothing parameter
      A_lambda = 1, B_lambda = 1,
      var_gamma = 100,
      A_tau2 = 3, B_tau2 = 2,  # IG(A,B) for tau^2 (shape A, scale B)
      A_sigma2 = 3, B_sigma2 = 2,  # IG(A,B) for sigma^2 (shape A, scale B)
      # MCMC
      n_iter = 20000,
      n_burn = 10000,
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
    surface_q025  <- apply(surface_draws, 2, quantile, 0.025, na.rm = TRUE)
    surface_q975  <- apply(surface_draws, 2, quantile, 0.975, na.rm = TRUE)
    
    ###############################################################
    ## surface estimation metrics
    ###############################################################
    
    beta_resid     <- surface_mean - beta_true_grid
    beta_ise       <- mean(beta_resid^2)           # approx integral over the 2D grid
    beta_bias_mean <- mean(beta_resid)             # signed mean bias
    beta_bias_max  <- max(abs(beta_resid))         # max absolute pointwise bias
    beta_cover_pw  <- mean(                        # pointwise 95% CrI coverage
      (beta_true_grid >= surface_q025) & (beta_true_grid <= surface_q975)
    )
    
    ###############################################################
    ## scalar parameter metrics
    ###############################################################
    
    # --- gamma (intercept + Z1 + Z2) ---
    # For Cox DGM: gamma enters on the log-hazard scale; the AFT model
    # estimates gamma on the log-time scale.  The two are related by
    # gamma_AFT ≈ -gamma_Cox / rho, so comparisons are approximate.
    gamma_true  <- sim_data$gamma
    gamma_est   <- colMeans(fit$gamma)
    gamma_bias  <- gamma_est  - gamma_true
    gamma_se    <- (gamma_est - gamma_true)^2
    q_gamma     <- apply(fit$gamma, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    gamma_cover <- (gamma_true >= q_gamma[1, ]) & (gamma_true <= q_gamma[2, ])
    
    # --- tau2 ---
    # Cox DGM: frailty u_j is on the log-hazard scale.
    # AFT model: frailty u_j is on the log-time scale.
    # Relationship for Weibull baseline: tau2_AFT_implied = tau2_Cox / rho^2.
    # Bias and coverage are computed against tau2_implied (the AFT-scale target).
    tau2_true    <- sim_data$tau^2
    tau2_implied <- if (as.character(family) == "cox") {
      tau2_true / sim_data$rho^2
    } else {
      tau2_true     # AFT DGM: same scale, no correction needed
    }
    tau2_est   <- mean(fit$tau2)
    tau2_bias  <- tau2_est - tau2_implied   # compare against implied target
    tau2_se    <- (tau2_est - tau2_implied)^2
    tau2_cover <- (tau2_implied >= quantile(fit$tau2, 0.025)) &
      (tau2_implied <= quantile(fit$tau2, 0.975))
    
    # --- sigma2 ---
    # For lognormal and Weibull AFT: sigma is the AFT error SD, so sigma2 = sigma^2
    # is directly estimated by the fitted model and meaningful to evaluate.
    # For Weibull specifically: sigma = 1/rho_weibull, so sigma2_true = 1/rho^2.
    # For Cox DGM: no AFT error SD exists -> NA.
    sigma2_true  <- if (as.character(family) == "cox") NA_real_ else sim_data$sigma^2
    sigma2_est   <- mean(fit$sigma2)
    sigma2_bias  <- sigma2_est - sigma2_true   # NA for Cox DGM
    sigma2_se    <- (sigma2_est - sigma2_true)^2
    sigma2_cover <- if (!is.na(sigma2_true)) {
      (sigma2_true >= quantile(fit$sigma2, 0.025)) &
        (sigma2_true <= quantile(fit$sigma2, 0.975))
    } else NA
    
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
      level          = 0.95
    )
    
    # Risk score: negative sum of predicted survival probabilities over tgrid.
    # Subjects with lower survival across all time points have a larger (less
    # negative) score, i.e. higher risk.  Consistent with reference code:
    #   risk_bayes <- -rowSums(S_bayes)
    risk_score <- -rowSums(pred_surv$mean)   # length N_test
    
    # Uno's IPCW C-index.
    # Upweights comparable pairs by the inverse probability of censoring,
    # making it consistent even under dependent censoring.
    # UnoC(Surv.rsp, Surv.rsp.new, lpnew) where lpnew = risk score for test set.
    Surv_train <- survival::Surv(data$Y, data$delta)  # training set (changes each iter)
    c_index    <- survAUC::UnoC(
      Surv.rsp     = Surv_train,
      Surv.rsp.new = Surv_test,
      lpnew        = risk_score
    )
    
    # IPCW Integrated Brier Score.
    # tgrid spans 5th-95th percentile of test EVENT times (precomputed outside
    # the loop) to include tails where lognormal and Weibull diverge most.
    ibs <- cal_IPCW_Brier(
      S_mat       = pred_surv$mean,
      time_test   = test_data$Y,
      event_test  = test_data$delta,
      time_train  = data$Y,
      event_train = data$delta,
      tgrid       = tgrid
    )
    
    ###############################################################
    ## assemble surface data frame  (n_grid rows, one per eval point)
    ###############################################################
    
    df_surface <- data.frame(
      scenario     = scenario,        # joins to result$info and result$coef
      iter         = iter,
      s            = grid_df$s_grid,
      age          = grid_df$age,
      age_scaled   = a_grid_scaled,
      surface_mean = surface_mean,
      surface_q025 = surface_q025,
      surface_q975 = surface_q975,
      beta_true    = beta_true_grid   # constant within scenario; stored per row
    )                                 # for self-contained downstream plotting
    
    ###############################################################
    ## assemble one-row result data frames
    ###############################################################
    
    # ---- coefficient / estimation summary ----
    df_coef <- data.frame(
      # surface
      beta_ise       = beta_ise,
      beta_bias_mean = beta_bias_mean,
      beta_bias_max  = beta_bias_max,
      beta_cover_pw  = beta_cover_pw,
      # frailty variance
      # tau2_true    : DGM-scale (log-hazard for Cox, log-time for AFT)
      # tau2_implied : AFT-equivalent target for bias/coverage
      #                = tau2_true        for AFT DGM
      #                = tau2_true / rho² for Cox DGM (scale correction)
      tau2_true    = tau2_true,
      tau2_implied = tau2_implied,
      tau2_est     = tau2_est,
      tau2_bias    = tau2_bias,       # tau2_est - tau2_implied
      tau2_se      = tau2_se,
      tau2_cover   = tau2_cover,      # I(tau2_implied in 95% CrI of fit$tau2)
      # residual variance (AFT only; NA for Cox DGM)
      sigma2_true  = sigma2_true,
      sigma2_est   = sigma2_est,
      sigma2_bias  = sigma2_bias,
      sigma2_se    = sigma2_se,
      sigma2_cover = sigma2_cover,
      # predictive — discrimination
      # C-index: P(pred_i > pred_j | T_i > T_j); > 0.5 = good model.
      # Insensitive to distributional misspecification (depends on LP ordering).
      c_index      = c_index,
      # predictive — calibration
      # IBS: integrates over 5th-95th percentile of test event times (tails included).
      # Sensitive to distributional misspecification; expect higher IBS for Cox DGM.
      ibs          = ibs
    )
    
    # gamma columns added dynamically (intercept = index 0, Z1 = 1, Z2 = 2)
    idx <- seq_along(gamma_est) - 1L
    df_coef[paste0("gamma_true_",  idx)] <- as.list(gamma_true)
    df_coef[paste0("gamma_est_",   idx)] <- as.list(gamma_est)
    df_coef[paste0("gamma_bias_",  idx)] <- as.list(gamma_bias)
    df_coef[paste0("gamma_se_",    idx)] <- as.list(gamma_se)
    df_coef[paste0("gamma_cover_", idx)] <- as.list(gamma_cover)
    
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
      tau         = tau,
      sigma       = if (as.character(family) == "cox") NA_real_ else sigma,
      rho         = if (as.character(family) == "weibull") 1/sigma
      else if (as.character(family) == "cox") rho
      else NA_real_,
      censor_rate = censor_rate,
      event_rate  = mean(data$delta),
      time        = time
    )
    
    list(info = df_info, coef = df_coef, surface = df_surface)
    
  }, error = function(e) {
    warning(sprintf(
      "Iteration %d skipped due to error:\n  %s",
      iter, e$message
    ))
    NULL
  })
  
  if (!is.null(res)) {
    info_list[[iter]]    <- res$info
    coef_list[[iter]]    <- res$coef
    surface_list[[iter]] <- res$surface
  }
} # end for loop

# drop NULL entries and bind into data frames
result <- list(
  info    = dplyr::bind_rows(Filter(Negate(is.null), info_list)),
  coef    = dplyr::bind_rows(Filter(Negate(is.null), coef_list)),
  surface = dplyr::bind_rows(Filter(Negate(is.null), surface_list))
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