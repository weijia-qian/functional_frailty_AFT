#' Calculate Integrated Brier Score (IPCW)
#' 
#' @param S_mat Matrix of predicted survival probabilities (Rows = test subjects, Cols = tgrid)
#' @param time_test Vector of observed survival times in the test set
#' @param event_test Vector of event indicators in the test set (1 = event, 0 = censored)
#' @param time_train Vector of observed survival times in the training set
#' @param event_train Vector of event indicators in the training set
#' @param tgrid Vector of time points matching the columns of S_mat
cal_Brier <- function(S_mat, time_test, event_test, time_train, event_train, tgrid) {
  
  N <- nrow(S_mat)
  
  # 1. Fit the Censoring Distribution G(t) = P(C > t) on the training data
  # Note: 1 - event_train means we are modeling the time UNTIL censoring
  G_fit <- survfit(Surv(time_train, 1 - event_train) ~ 1)
  
  # Helper function to safely extract G(t) using R's built-in step function logic
  get_G <- function(t_eval) {
    # 'extend = TRUE' safely carries the last observation forward 
    # if test times exceed training times
    res <- summary(G_fit, times = t_eval, extend = TRUE)$surv
    
    # Prevent division by zero if the Kaplan-Meier curve hits exactly 0
    res[res == 0] <- 1e-5 
    return(res)
  }
  
  brier_scores <- numeric(length(tgrid))
  
  # 2. Calculate IPCW Brier Score at each time point in tgrid
  for (k in seq_along(tgrid)) {
    t_star <- tgrid[k]
    pred_surv <- S_mat[, k]
    
    # Condition 1: Subject had an event BEFORE OR AT t_star
    cond1 <- (time_test <= t_star) & (event_test == 1)
    
    # Condition 2: Subject survived PAST t_star
    cond2 <- (time_test > t_star)
    
    # Note: Subjects censored before t_star fall into neither condition 
    # and contribute 0 to the score (they are factored out by the IPCW weights of others).
    
    # Get censoring probabilities (weights)
    G_Ti <- get_G(time_test) # G evaluated at the subject's actual event time
    G_t_star <- get_G(t_star) # G evaluated at the prediction time point
    
    # Brier score components
    term1 <- cond1 * ((0 - pred_surv)^2) / G_Ti
    term2 <- cond2 * ((1 - pred_surv)^2) / G_t_star
    
    # Sum and average
    brier_scores[k] <- sum(term1 + term2, na.rm = TRUE) / N
  }
  
  # 3. Integrate over tgrid using the Trapezoidal rule
  dt <- diff(tgrid)
  ibs <- sum(dt * (brier_scores[-1] + brier_scores[-length(brier_scores)]) / 2)
  
  # Normalize by the maximum time to get the average Integrated Brier Score
  ibs <- ibs / max(tgrid)
  
  return(ibs)
}

cal_IPCW_Brier <- function(S_mat, time_test, event_test, time_train, event_train, tgrid) {
  n_test <- length(time_test)
  n_grid <- length(tgrid)
  
  # 1. Estimate Censoring Distribution G(t) using Kaplan-Meier on TRAINING data
  # Note: Reverse the event indicator (1 = censored, 0 = event)
  km_cens <- survival::survfit(survival::Surv(time_train, 1 - event_train) ~ 1)
  
  # Function to get G(t) safely using a step function
  get_G <- function(t_vals) {
    sf <- stepfun(km_cens$time, c(1, km_cens$surv))
    G_vals <- sf(t_vals)
    # Prevent division by zero if evaluation time exceeds last observed censoring
    G_vals[G_vals == 0] <- min(km_cens$surv[km_cens$surv > 0]) 
    return(G_vals)
  }
  
  BS_t <- numeric(n_grid)
  
  for (k in 1:n_grid) {
    t_star <- tgrid[k]
    pred_S <- S_mat[, k]
    
    # Term 1: Subject had event before or at t_star (and wasn't censored before event)
    # Weight uses G(T_i)
    term1_idx <- (time_test <= t_star) & (event_test == 1)
    G_Ti <- get_G(time_test[term1_idx])
    term1 <- (pred_S[term1_idx]^2) / G_Ti
    
    # Term 2: Subject survived past t_star
    # Weight uses G(t_star)
    term2_idx <- (time_test > t_star)
    G_tstar <- get_G(t_star)
    term2 <- ((1 - pred_S[term2_idx])^2) / G_tstar
    
    # Brier score at time t_star
    BS_t[k] <- (sum(term1) + sum(term2)) / n_test
  }
  
  # Integrated Brier Score using trapezoidal rule
  dt <- diff(tgrid)
  IBS <- sum(dt * (BS_t[-1] + BS_t[-length(BS_t)]) / 2) / (max(tgrid) - min(tgrid))
  
  return(IBS)
}
