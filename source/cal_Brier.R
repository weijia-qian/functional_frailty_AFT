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

#' Calculate the IPCW Integrated Brier Score
#'
#' @param S_mat       Predicted survival probabilities, rows = test subjects,
#'                    columns = tgrid.
#' @param time_test,event_test    Observed time / event indicator, test set.
#' @param time_train,event_train  Observed time / event indicator, training set.
#'                    Used only when `G_data = "train"`.
#' @param tgrid       Time points matching the columns of `S_mat`.
#' @param G_data      Which sample to estimate the censoring distribution G from.
#'                    "test" (default) follows Graf et al. (1999) / Gerds &
#'                    Schumacher (2006): the IPCW correction applies to the
#'                    observations being weighted, so G belongs to the evaluation
#'                    sample.  "train" reproduces the previous behaviour.
#' @param eps_G       Lower bound on G.  Grid points where G(t) < eps_G are
#'                    dropped (the IPCW weight is not estimable there) and the
#'                    remaining weights are clamped at eps_G as a guard.
#'
#' @details
#' Estimating G on the training data while `tgrid` is built from test event
#' times was the source of extreme IBS values.  Because Y = min(T, C) with a
#' bounded censoring time, the largest training observation is essentially
#' always a censoring time, so the reverse Kaplan-Meier hits exactly 0 at
#' max(time_train).  Whenever the upper end of tgrid reached past that point,
#' G(t) = 0 and the old floor -- `G[G == 0] <- min(surv[surv > 0])`, i.e. the
#' smallest positive KM value -- produced weights of 1/tiny and IBS values up to
#' ~57.  Scoring oracle predictions, blow-ups occurred in 2-3 of 40 replicates
#' at 75% censoring in every family, and IBS > 1 exactly when G hit 0.
#'
#' Taking G from the test set removes the mismatch (tgrid is built from the test
#' sample, so it is inside its own support by construction) and, because the
#' test set is fixed across the iterations of a scenario, it also holds the
#' weights and the evaluation horizon constant so the replicates estimate the
#' same functional.  The event term now uses the left limit G(T_i-), which is
#' the standard IPCW weight.
#'
#' @return The integrated Brier score, with attributes `horizon` (the largest
#'   retained grid point), `n_kept` and `BS_t` (the pointwise scores).
cal_IPCW_Brier <- function(S_mat, time_test, event_test, time_train, event_train,
                           tgrid, G_data = c("test", "train"), eps_G = 0.05) {

  G_data  <- match.arg(G_data)
  g_time  <- if (G_data == "test") time_test  else time_train
  g_event <- if (G_data == "test") event_test else event_train

  # 1. Censoring distribution G(t) = P(C > t); reverse the event indicator.
  km_cens <- survival::survfit(survival::Surv(g_time, 1 - g_event) ~ 1)
  G_right <- stats::stepfun(km_cens$time, c(1, km_cens$surv))                # G(t)
  G_left  <- stats::stepfun(km_cens$time, c(1, km_cens$surv), right = TRUE)  # G(t-)

  # 2. Restrict to the range where the IPCW weight is actually estimable.
  keep <- G_right(tgrid) >= eps_G
  if (sum(keep) < 2L) return(NA_real_)
  tgrid <- tgrid[keep]
  S_mat <- S_mat[, keep, drop = FALSE]

  n_test <- length(time_test)
  BS_t   <- numeric(length(tgrid))

  # 3. IPCW Brier score at each retained time point.
  for (k in seq_along(tgrid)) {
    t_star <- tgrid[k]
    pred_S <- S_mat[, k]

    # Term 1: event at or before t_star, weighted by G(T_i-)
    idx1 <- (time_test <= t_star) & (event_test == 1)
    w1   <- pmax(G_left(time_test[idx1]), eps_G)
    term1 <- (pred_S[idx1]^2) / w1

    # Term 2: still at risk past t_star, weighted by G(t_star)
    idx2 <- time_test > t_star
    w2   <- pmax(G_right(t_star), eps_G)
    term2 <- ((1 - pred_S[idx2])^2) / w2

    # Subjects censored before t_star fall into neither term.
    BS_t[k] <- (sum(term1) + sum(term2)) / n_test
  }

  # 4. Integrate over the retained grid (trapezoidal), normalised by its span.
  dt  <- diff(tgrid)
  IBS <- sum(dt * (BS_t[-1] + BS_t[-length(BS_t)]) / 2) / (max(tgrid) - min(tgrid))

  structure(IBS, horizon = max(tgrid), n_kept = length(tgrid), BS_t = BS_t)
}
