####################################################################
# Weijia Qian
# Jan 10, 2026
#
# This file simulates survival data under different data generation mechanisms
# and fits functional frailty AFT models
####################################################################

suppressPackageStartupMessages(library(here))
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
source(here("source", "simulate_AFT.R"))
load(here("data", "dat_func.Rdata")) # real data

###############################################################
## set simulation design elements
###############################################################
family = c("lognormal")
n_cluster = c(100, 200, 500)
n_subject = c(5)
nS = c(100)
beta_type = c('monotone', 'peak1', 'peak2', 'wavy')
tau = c(0.5, 1, 2)
sigma = c(0.5, 1, 2)
censor_rate = c(0.25)
N_iter = 500

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
  scenario = 1
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
  sim_data <- simulate_AFT(family = as.character(family),  # "lognormal" or "loglogistic"
                           n_cluster = n_cluster,       # number of clusters
                           n_subject = n_subject,        # number of subjects per cluster
                           nS = nS,              # density of s grid
                           k0 = 20,
                           alpha = 0,
                           beta_type = as.character(beta_type),    # "monotone"/"peak1"/"peak2"/"wavy"
                           gamma = c(0.5, 0.3, -0.2),  # c(intercept, gamma1, gamma2)
                           tau = tau,               # SD of random noise
                           sigma = sigma,             # SD of frailty term
                           censor_rate = censor_rate      
                           )
  
  res <- tryCatch({
    
    ###############################################################
    ## fit functional frailty AFT  model
    ###############################################################
    
    # extract elements from simulated data
    data <- sim_data$data
    s_grid = sim_data$coefficients$time
    Z <- model.matrix(~ Z1 + Z2, data = data)
    
    # run Gibbs sampler
    tic()
    fit <- gibbs_functional_frailty(
      time = data$Y,
      status = data$delta,
      cluster_id = data$cluster_id,
      Z = Z,
      X = data$X,
      s_grid = s_grid,
      # tuning / priors
      K = 20,
      a_pen = 0.001,      # MUST be > 0 to make D PD
      lambda_init = 1000,         # smoothing parameter
      A_lambda = 1, B_lambda = 0.001,
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
    ## integrated mean squared errors, pointwise CIs
    ###############################################################
    beta_true <- sim_data$coefficients$beta
    beta_est <- fit$beta_mean
    beta_bias <- beta_est - beta_true
    beta_ise <- mean((beta_est - beta_true)^2)
    beta_cover <- mean((beta_true >= fit$beta_q025) & (beta_true <= fit$beta_q975))
    
    gamma_true <- sim_data$coefficients$gamma[1,]
    gamma_est <- colMeans(fit$gamma)
    gamma_bias <- gamma_est - gamma_true
    gamma_se <- (gamma_est - gamma_true)^2
    q <- apply(fit$gamma, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    gamma_cover <- (gamma_true >= q[1, ]) & (gamma_true <= q[2, ])
    
    tau2_true <- (sim_data$coefficients$tau[1])^2
    tau2_est <- mean(fit$tau2)
    tau2_bias <- tau2_est - tau2_true
    tau2_se <- (tau2_est - tau2_true)^2
    tau2_cover <- (tau2_true >= quantile(fit$tau2, 0.025)) & (tau2_true <= quantile(fit$tau2, 0.975))
    
    sigma2_true <- (sim_data$coefficients$sigma[1])^2
    sigma2_est <- mean(fit$sigma2)
    sigma2_bias <- sigma2_est - sigma2_true
    sigma2_se <- (sigma2_est - sigma2_true)^2
    sigma2_cover <- (sigma2_true >= quantile(fit$sigma2, 0.025)) & (sigma2_true <= quantile(fit$sigma2, 0.975))
    
    # ---- one-row coef summary ----
    df_coef <- data.frame(
      beta_ise   = beta_ise,
      beta_cover = beta_cover,
      tau2_est   = tau2_est,
      tau2_bias  = tau2_bias,
      tau2_se    = tau2_se,
      tau2_cover = tau2_cover,
      sigma2_est   = sigma2_est,
      sigma2_bias  = sigma2_bias,
      sigma2_se    = sigma2_se,
      sigma2_cover = sigma2_cover
    )
    
    df_coef[paste0("gamma_est_",   seq_along(gamma_est) - 1)]    <- as.list(gamma_est)
    df_coef[paste0("gamma_bias_",  seq_along(gamma_bias) - 1)]   <- as.list(gamma_bias)
    df_coef[paste0("gamma_se_",    seq_along(gamma_se) - 1)]     <- as.list(gamma_se)
    df_coef[paste0("gamma_cover_", seq_along(gamma_cover) - 1)]  <- as.list(gamma_cover)
    df_coef[paste0("beta_bias_", seq_along(beta_bias))]  <- as.list(beta_bias)
    
    # ---- one-row info ----
    df_info <- data.frame(
      scenario    = scenario,
      iter        = iter,
      seed        = seeds[[iter]],
      family      = family,
      n_cluster   = n_cluster,
      n_subject   = n_subject,
      nS          = nS,
      beta_type   = beta_type,
      tau         = tau,
      sigma       = sigma,
      censor_rate = censor_rate,
      event_rate  = mean(data$delta),
      time        = time
    ) 
    
    # ---- one-row beta bias summary ----
    
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
