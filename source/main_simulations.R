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
family = c("lognormal", "loglogistic")
n_cluster = c(50, 100, 200)
n_subject = c(5)
nS = c(100)
beta_type = c('simple')
tau = c(0.1, 0.5, 2)
seed_start = 1000
N_iter = 100

params = expand.grid(family = family,
                     n_cluster = n_cluster,
                     n_subject = n_subject,
                     nS = nS,
                     beta_type = beta_type,
                     tau = tau)

## define number of simulations and parameter scenarios
if(doLocal) {
  scenario = 10
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

###############################################################
## run simulations
###############################################################
results = vector("list", length = N_iter)

for(iter in 1:N_iter){
  print(iter)
  # set seed
  seed.iter = (seed_start - 1) * N_iter + iter

  # simulate data
  sim_data <- simulate_AFT(family = as.character(family), n_cluster = n_cluster, n_subject = n_subject,
                             nS = nS, beta_type = as.character(beta_type), tau = tau, seed = seed.iter)

  res <- tryCatch({
  ###############################################################
  ## fit functional frailty AFT  model
  ###############################################################
  tic()
    
  # extract elements from simulated data
  data <- sim_data$data
  s_grid = sim_data$coefficients$time
  Z <- model.matrix(~ Z1 + Z2, data = data)
    
  # run Gibbs sampler
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
      lambda = 5000,      
      sigma_gamma2 = 25, 
      A = 2, B = 1,      
      # MCMC
      n_iter = 10000,
      n_burn = 5000,
      n_thin = 1,
      seed = 42,
      verbose = TRUE
  )  
  
  time_stamp <- toc(quiet = TRUE)
  time <- time_stamp$toc - time_stamp$tic
  
  ###############################################################
  ## integrated mean squared errors, pointwise CIs
  ###############################################################
  # calculate MSE/ISME and coverage
  beta1_true <- sim_data$coefficients$beta1
  beta1_imse <- mean((beta1_true - fit$beta_mean)^2)
  beta1_cover <- mean((beta1_true > fit$beta_q025) & (beta1_true < fit$beta_q975))
  
  tau2_true <- (sim_data$coefficients$tau[1])^2
  tau2_mse <- (tau2_true - quantile(fit$tau2, 0.5))^2
  tau2_cover <- (tau2_true > quantile(fit$tau2, 0.025)) & (tau2_true < quantile(fit$tau2, 0.975))
  
  sigma2_true <- (sim_data$coefficients$sigma[1])^2
  sigma2_mse <- (sigma2_true - quantile(fit$sigma2, 0.5))^2
  sigma2_cover <- (sigma2_true > quantile(fit$sigma2, 0.025)) & (sigma2_true < quantile(fit$sigma2, 0.975))
  
  # summary
  df_coef <- data.frame(beta1_imse, beta1_cover,
                        tau2_mse, tau2_cover,
                        sigma2_mse, sigma2_cover)
  
  df_info <- data.frame(scenario = scenario,
                        iter = iter,
                        seed = seed.iter,
                        family = family,
                        n_cluster = n_cluster,
                        n_subject = n_subject,
                        nS = nS,
                        beta_type = beta_type,
                        tau = tau,
                        censor_rate = 1 - mean(data$delta),
                        time
                        )
  
  list(
    info = df_info,
    coef = df_coef
  )
  
  }, error = function(e) {
    warning(sprintf(
      "Iteration %d skipped due to error:\n  %s",
      iter, e$message
    ))
    NULL
  })
  
  ## only save non-NULL results
  if (!is.null(res)) {
    results[[iter]] <- res
  }

} # end for loop

# drop NULL entries
results <- Filter(Negate(is.null), results)

# record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())
dir.create(file.path(here::here("outputs"), Date), showWarnings = FALSE)

filename = paste0(here::here("outputs", Date), "/", scenario, ".RDA")
save(results,
     file = filename)

###############################################################
## end sim
###############################################################


