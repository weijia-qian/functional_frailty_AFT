####################################################################
# Weijia Qian
# Jan 10, 2026
#
# Simulate survival data under a functional frailty AFT model with a
# bivariate coefficient surface beta(s, age), then fit
# gibbs_frailty_interaction() and extract estimation metrics.
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
#   tau2_{est,bias,se,cover}
#   sigma2_{est,bias,se,cover}
#
# Computational:
#   time             wall-clock seconds for the Gibbs sampler
####################################################################

suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(refund))
suppressPackageStartupMessages(library(splines))
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
source(here("source", "simulate_AFT_interaction.R"))

###############################################################
## set simulation design elements
###############################################################
family = c("lognormal")
n_cluster = c(25, 100, 300)
n_subject = c(5)
# n_subject = c(5, 10, 50)
nS = c(100)
beta_type = c("additive", "bilinear", "ridge")
tau = c(0.1, 0.5, 1)
sigma = c(0.5, 1, 2)
censor_rate = c(0.1, 0.25, 0.5)
N_iter = 100

params = expand.grid(family = family,
                     n_cluster = n_cluster,
                     n_subject = n_subject,
                     nS = nS,
                     beta_type = beta_type,
                     tau = tau,
                     sigma = sigma,
                     censor_rate = censor_rate)

## define number of simulations and parameter scenarios
if(doLocal) {
  scenario = 293
  N_iter = 2
}else{
  # defined from batch script params
  scenario <- as.numeric(commandArgs(trailingOnly=TRUE))
}

family = params$family[scenario]
n_cluster = params$n_cluster[scenario]
n_subject = params$n_subject[scenario]
nS = params$nS[scenario]
beta_type = params$beta_type[scenario]
tau = params$tau[scenario]
sigma = params$sigma[scenario]
censor_rate = params$censor_rate[scenario]

###############################################################
## run simulations
###############################################################

# collect one-row data frames per successful iteration
coef_list <- vector("list", length = N_iter)
info_list <- vector("list", length = N_iter)
seeds <- sample.int(1e8, N_iter)

for(iter in 1:N_iter){
  cat("Scenario", scenario, "| Iteration", iter, "\n")
  # set seed
  set.seed(seeds[[iter]])
  
  # simulate data
  sim_data <- simulate_AFT_interaction(
    family       = as.character(family),
    n_cluster    = n_cluster,
    n_subject    = n_subject,
    nS           = nS,
    k0           = 20,
    alpha        = 0.7,
    beta_type = as.character(beta_type),   # "additive" or "bilinear"
    gamma        = c(0.5, 0.3, -0.2),
    tau          = tau,
    sigma        = sigma,
    censor_rate  = censor_rate
  )
  
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
      K_s = 10,
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
    ## reconstruct posterior surface beta(s, a) on an eval grid
    ## fit$b holds S x K_total draws; the surface is b %*% t(Phi)
    ## where Phi = PredictMat(sm_tensor, grid_df)
    ###############################################################
    
    # Evaluation grid: 30 s-points x 20 age-points
    # Use raw age (same scale as training) for PredictMat;
    # use scaled age a in [0,1] for the true surface beta_fn(s, a).
    n_eval_s   <- 30
    n_eval_age <- 30
    age_range  <- sim_data$age_range
    s_eval     <- seq(min(s_grid), max(s_grid), length.out = n_eval_s)
    age_eval   <- seq(age_range[1], age_range[2], length.out = n_eval_age)
    a_eval     <- (age_eval - age_range[1]) / diff(age_range)  # scaled to [0,1]
    
    grid_df    <- expand.grid(s_grid = s_eval, age = age_eval)  # (n_eval_s * n_eval_age) rows
    
    # Basis matrix on the eval grid (reusing the sm_tensor stored in the fit)
    Phi_grid   <- PredictMat(fit$meta$sm_tensor, grid_df)       # n_grid x K_total
    
    # Posterior draws of the surface: S x n_grid
    surface_draws <- fit$b %*% t(Phi_grid)
    
    surface_mean  <- colMeans(surface_draws)
    surface_q025  <- apply(surface_draws, 2, quantile, 0.025, na.rm = TRUE)
    surface_q975  <- apply(surface_draws, 2, quantile, 0.975, na.rm = TRUE)
    
    # True surface values on the same grid
    a_grid_scaled <- (grid_df$age - age_range[1]) / diff(age_range)
    beta_true_grid <- mapply(sim_data$beta_fn, grid_df$s_grid, a_grid_scaled)
    
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
    gamma_true  <- sim_data$gamma                              # length-3 vector from DGM
    gamma_est   <- colMeans(fit$gamma)
    gamma_bias  <- gamma_est  - gamma_true
    gamma_se    <- (gamma_est - gamma_true)^2
    q_gamma     <- apply(fit$gamma, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    gamma_cover <- (gamma_true >= q_gamma[1, ]) & (gamma_true <= q_gamma[2, ])
    
    # --- tau2 ---
    tau2_true  <- sim_data$tau^2
    tau2_est   <- mean(fit$tau2)
    tau2_bias  <- tau2_est - tau2_true
    tau2_se    <- (tau2_est - tau2_true)^2
    tau2_cover <- (tau2_true >= quantile(fit$tau2, 0.025)) &
      (tau2_true <= quantile(fit$tau2, 0.975))
    
    # --- sigma2 ---
    sigma2_true  <- sim_data$sigma^2
    sigma2_est   <- mean(fit$sigma2)
    sigma2_bias  <- sigma2_est - sigma2_true
    sigma2_se    <- (sigma2_est - sigma2_true)^2
    sigma2_cover <- (sigma2_true >= quantile(fit$sigma2, 0.025)) &
      (sigma2_true <= quantile(fit$sigma2, 0.975))
    
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
      # variance components
      tau2_true    = tau2_true,
      tau2_est     = tau2_est,
      tau2_bias    = tau2_bias,
      tau2_se      = tau2_se,
      tau2_cover   = tau2_cover,
      sigma2_true  = sigma2_true,
      sigma2_est   = sigma2_est,
      sigma2_bias  = sigma2_bias,
      sigma2_se    = sigma2_se,
      sigma2_cover = sigma2_cover
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
      n_cluster   = n_cluster,
      n_subject   = n_subject,
      nS          = nS,
      beta_type = as.character(beta_type),
      tau         = tau,
      sigma       = sigma,
      censor_rate = censor_rate,
      event_rate  = mean(data$delta),
      time        = time
    )
    
    list(info = df_info, coef = df_coef)
    
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
  }
} # end for loop

# drop NULL entries and bind into data frames
result <- list(
  info = dplyr::bind_rows(Filter(Negate(is.null), info_list)),
  coef = dplyr::bind_rows(Filter(Negate(is.null), coef_list))
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
