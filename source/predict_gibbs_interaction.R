#' Predict from a Bayesian Functional Interaction AFT Model with Frailty
predict_gibbs_interaction <- function(fit,
                                      X_new,
                                      age_new,
                                      Z_new = NULL,
                                      cluster_id_new = NULL,
                                      type = c("link", "response", "survival"),
                                      times = NULL,
                                      level = 0.95) {
  
  type <- match.arg(type)
  require(mgcv)
  
  # ---- 1. Basic Checks and Meta Extraction ----
  if (is.null(fit$b) || is.null(fit$gamma)) {
    stop("`fit` does not appear to be a valid gibbs_frailty_interaction output.")
  }
  
  meta <- fit$meta
  s_grid <- meta$s_grid
  sm_tensor <- meta$sm_tensor
  
  if (is.null(s_grid)) stop("`fit$meta$s_grid` is missing.")
  if (is.null(sm_tensor)) stop("`fit$meta$sm_tensor` is missing. Ensure the model saved it.")
  
  N_new <- nrow(X_new)
  nS_fit <- length(s_grid)
  nS_new <- ncol(X_new)
  
  if (nS_fit != nS_new) {
    stop(sprintf("Number of grid points in fit (%d) does not match ncol(X_new) = %d.", 
                 nS_fit, nS_new))
  }
  if (length(age_new) != N_new) {
    stop("length(age_new) must equal nrow(X_new).")
  }
  
  S_iter <- nrow(fit$b) 
  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)
  
  # ---- 2. Integration Weights ----
  P <- length(s_grid)
  ds <- diff(s_grid)
  w_dt <- numeric(P)
  w_dt[1] <- ds[1] / 2
  w_dt[P] <- ds[P - 1] / 2
  if (P > 2) w_dt[2:(P - 1)] <- (ds[1:(P - 2)] + ds[2:(P - 1)]) / 2
  
  # ---- 3. Calculate Integrated Interaction Draws (N_new x S_iter) ----
  K_total <- meta$K_total
  W_new <- matrix(0, nrow = N_new, ncol = K_total)
  
  # Construct the new integrated tensor matrix W_new
  for (i in seq_len(N_new)) {
    # Evaluate the tensor basis at the functional grid AND the new subject's age
    df_i <- data.frame(s_grid = s_grid, age = rep(age_new[i], P))
    Phi_i <- PredictMat(sm_tensor, df_i) # Size: P x K_total
    
    # Integrate: sum_t( X_new_i(t) * Phi_i(t, age_new_i) * w_dt )
    W_new[i, ] <- colSums(X_new[i, ] * Phi_i * w_dt)
  }
  
  # Calculate draws for the interaction term
  int_draws <- W_new %*% t(fit$b) 
  
  # ---- 4. Calculate Scalar Component Draws (N_new x S_iter) ----
  if (!is.null(Z_new) && meta$M > 0) {
    if (ncol(Z_new) != meta$M) stop("Z_new dimensions do not match training data.")
    scalar_draws <- Z_new %*% t(fit$gamma)
  } else {
    scalar_draws <- matrix(0, nrow = N_new, ncol = S_iter)
  }
  
  # ---- 5. Calculate Frailty Component Draws (N_new x S_iter) ----
  frailty_draws <- matrix(0, nrow = N_new, ncol = S_iter)
  if (!is.null(cluster_id_new) && meta$J > 0) {
    
    # Safely map categorical factors / character strings to integer column indices
    if (!is.null(meta$cluster_levels)) {
      cid_indices <- match(as.character(cluster_id_new), meta$cluster_levels)
    } else {
      cid_indices <- suppressWarnings(as.numeric(as.character(cluster_id_new)))
      if (any(is.na(cid_indices) & !is.na(cluster_id_new))) {
        if (is.factor(cluster_id_new)) {
          cid_indices <- as.integer(cluster_id_new)
        }
      }
    }
    
    # Loop over the new numeric indices
    for (i in seq_len(N_new)) {
      idx <- cid_indices[i]
      if (!is.na(idx) && idx > 0 && idx <= meta$J) {
        frailty_draws[i, ] <- fit$u[, idx]
      }
    }
  }
  
  # ---- 6. Combine for Linear Predictor (mu) Draws ----
  mu_draws <- int_draws + scalar_draws + frailty_draws
  
  # ---- 7. Process Output by Type ----
  res <- list()
  
  if (type == "link") {
    res$mean  <- rowMeans(mu_draws)
    res$lower <- apply(mu_draws, 1, quantile, probs[1])
    res$upper <- apply(mu_draws, 1, quantile, probs[2])
    
  } else if (type == "response") {
    resp_draws <- exp(mu_draws)
    res$mean  <- rowMeans(resp_draws)
    res$lower <- apply(resp_draws, 1, quantile, probs[1])
    res$upper <- apply(resp_draws, 1, quantile, probs[2])
    
  } else { # type == "survival"
    if (is.null(times) || any(times <= 0)) {
      stop("`times` must be a numeric vector of strictly positive values.")
    }
    
    sigma_draws <- sqrt(fit$sigma2) 
    
    S_mean  <- matrix(NA, nrow = N_new, ncol = length(times))
    S_lower <- matrix(NA, nrow = N_new, ncol = length(times))
    S_upper <- matrix(NA, nrow = N_new, ncol = length(times))
    
    for (k in seq_along(times)) {
      lt <- log(times[k])
      z_num <- lt - mu_draws 
      z_draws <- sweep(z_num, 2, sigma_draws, `/`)
      
      S_draws <- 1 - pnorm(z_draws)
      
      S_mean[, k]  <- rowMeans(S_draws)
      S_lower[, k] <- apply(S_draws, 1, quantile, probs[1])
      S_upper[, k] <- apply(S_draws, 1, quantile, probs[2])
    }
    
    colnames(S_mean) <- colnames(S_lower) <- colnames(S_upper) <- as.character(times)
    
    res$mean  <- S_mean
    res$lower <- S_lower
    res$upper <- S_upper
  }
  
  if (!is.null(rownames(X_new))) {
    if (type == "survival") {
      rownames(res$mean) <- rownames(res$lower) <- rownames(res$upper) <- rownames(X_new)
    } else {
      names(res$mean) <- names(res$lower) <- names(res$upper) <- rownames(X_new)
    }
  }
  
  return(res)
}