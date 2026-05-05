#' Predict from a Bayesian Linear Functional AFT Model with Frailty (with CrI)
predict_gibbs_frailty <- function(fit,
                                  X_new,
                                  Z_new = NULL,
                                  cluster_id_new = NULL,
                                  type = c("link", "response", "survival"),
                                  times = NULL,
                                  level = 0.95) {
  
  type <- match.arg(type)
  
  # ---- 1. Basic Checks and Meta Extraction ----
  if (is.null(fit$b) || is.null(fit$gamma)) {
    stop("`fit` does not appear to be a valid gibbs_functional_frailty output.")
  }
  
  meta <- fit$meta
  s_grid <- meta$s_grid
  if (is.null(s_grid)) stop("`fit$meta$s_grid` is missing.")
  
  N_new <- nrow(X_new)
  nS_fit <- length(s_grid)
  nS_new <- ncol(X_new)
  
  if (nS_fit != nS_new) {
    stop(sprintf("Number of grid points in fit (%d) does not match ncol(X_new) = %d.", 
                 nS_fit, nS_new))
  }
  
  S_iter <- nrow(fit$b) 
  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)
  
  # ---- 2. FIXED: Rebuild the Basis Matrix (Bmat) accurately ----
  K_fit <- ncol(fit$b) # Exact number of columns we need
  
  if (meta$basis_type == "bs") {
    Bmat <- splines::bs(s_grid, df = K_fit, intercept = TRUE)  
  } else {
    require(mgcv)
    # mgcv dynamically drops dimensions for certain constraints.
    # We loop to find the exact 'k' that recreates the training dimensions.
    Bmat <- NULL
    for (k_try in K_fit:(K_fit + 5)) {
      sm_obj <- smoothCon(s(s_grid, k = k_try, bs = meta$basis_type), 
                          data = data.frame(s_grid = s_grid), 
                          absorb.cons = FALSE)[[1]]
      if (ncol(sm_obj$X) == K_fit) {
        Bmat <- sm_obj$X
        break
      }
    }
    if (is.null(Bmat)) {
      stop("Could not reconstruct the mgcv basis matrix to match training dimensions.")
    }
  }
  
  # ---- 3. Calculate Functional Component Draws (N_new x S_iter) ----
  P <- length(s_grid)
  ds <- diff(s_grid)
  w <- numeric(P)
  w[1] <- ds[1] / 2
  w[P] <- ds[P - 1] / 2
  if (P > 2) w[2:(P - 1)] <- (ds[1:(P - 2)] + ds[2:(P - 1)]) / 2
  
  Xw_new <- sweep(X_new, 2, w, `*`)
  
  # Now conformable!
  beta_draws <- fit$b %*% t(Bmat) 
  func_draws <- Xw_new %*% t(beta_draws)
  
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
    
    # NEW: Safely map categorical factors / character strings to integer column indices
    if (!is.null(meta$cluster_levels)) {
      # Best case scenario: The training function saved the original cluster levels
      cid_indices <- match(as.character(cluster_id_new), meta$cluster_levels)
    } else {
      # Fallback: Extract numeric index from strings like "1", "2"
      cid_indices <- suppressWarnings(as.numeric(as.character(cluster_id_new)))
      
      # If they are arbitrary factors ("A", "B") and couldn't be parsed to numerics,
      # fall back to the underlying integer mapping of the factor.
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
  mu_draws <- func_draws + scalar_draws + frailty_draws
  
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