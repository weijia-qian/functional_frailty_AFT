#' Predict from a Bayesian Linear Functional AFT Model with Frailty (with CrI)
predict_gibbs_frailty <- function(fit,
                                  X_new,
                                  Z_new = NULL,
                                  cluster_id_new = NULL,
                                  type = c("link", "response", "survival"),
                                  times = NULL,
                                  level = 0.95,
                                  quantiles = TRUE) {
  # `quantiles = FALSE` returns only $mean and skips the pointwise credible
  # bands, which dominate the runtime at large n_draws.  Mirrors the same
  # argument on predict_gibbs_interaction().
  
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
  
  # beta(s) draws -> integrated functional design, N_new x K
  W_new <- Xw_new %*% Bmat

  # ---- 4-7. Scalar + frailty components, combined and summarised in ROW BLOCKS ----
  # int/scalar/frailty draws are each N_new x S_iter (76 MB at N_new = 1000 and
  # 10,000 draws), and the survival branch forms three more of that size on every
  # time step -- ~1.2 GB peak.  Processing subjects in blocks caps the transient
  # at block x S_iter without changing any result: every output is a per-subject
  # summary, so the blocking only changes the order of allocation.
  block <- max(1L, min(N_new, floor(2e6 / max(S_iter, 1L))))

  cid_indices <- NULL
  if (!is.null(cluster_id_new) && meta$J > 0) {
    if (!is.null(meta$cluster_levels)) {
      cid_indices <- match(as.character(cluster_id_new), meta$cluster_levels)
    } else {
      cid_indices <- suppressWarnings(as.numeric(as.character(cluster_id_new)))
      if (any(is.na(cid_indices) & !is.na(cluster_id_new)) && is.factor(cluster_id_new)) {
        cid_indices <- as.integer(cluster_id_new)
      }
    }
  }
  if (!is.null(Z_new) && meta$M > 0 && ncol(Z_new) != meta$M)
    stop("Z_new dimensions do not match training data.")

  # Transposes hoisted out of the block loop: t(fit$b) is K x S_iter (11 MB at
  # 10,000 draws) and was otherwise rebuilt for every block.
  bt <- t(fit$b)
  gt <- if (!is.null(Z_new) && meta$M > 0) t(fit$gamma) else NULL

  # mu draws for one block of subjects
  mu_block <- function(idx) {
    mu <- W_new[idx, , drop = FALSE] %*% bt
    if (!is.null(gt)) mu <- mu + Z_new[idx, , drop = FALSE] %*% gt
    if (!is.null(cid_indices)) {
      j <- cid_indices[idx]
      ok <- !is.na(j) & j > 0 & j <= meta$J
      if (any(ok)) mu[ok, ] <- mu[ok, , drop = FALSE] + t(fit$u[, j[ok], drop = FALSE])
    }
    mu
  }

  res <- list()
  if (type %in% c("link", "response")) {
    res$mean <- numeric(N_new)
    if (quantiles) { res$lower <- numeric(N_new); res$upper <- numeric(N_new) }
    for (st in seq(1L, N_new, by = block)) {
      idx <- st:min(st + block - 1L, N_new)
      d <- mu_block(idx)
      if (type == "response") d <- exp(d)
      res$mean[idx] <- rowMeans(d)
      if (quantiles) {
        res$lower[idx] <- apply(d, 1, quantile, probs[1])
        res$upper[idx] <- apply(d, 1, quantile, probs[2])
      }
    }
  } else { # survival
    if (is.null(times) || any(times <= 0))
      stop("`times` must be a numeric vector of strictly positive values.")
    sigma_draws <- sqrt(fit$sigma2)
    nt <- length(times)
    res$mean <- matrix(NA_real_, N_new, nt)
    if (quantiles) {
      res$lower <- matrix(NA_real_, N_new, nt); res$upper <- matrix(NA_real_, N_new, nt)
    }
    for (st in seq(1L, N_new, by = block)) {
      idx <- st:min(st + block - 1L, N_new)
      mu <- mu_block(idx)
      for (k in seq_len(nt)) {
        S_draws <- 1 - pnorm(sweep(log(times[k]) - mu, 2, sigma_draws, `/`))
        res$mean[idx, k] <- rowMeans(S_draws)
        if (quantiles) {
          res$lower[idx, k] <- apply(S_draws, 1, quantile, probs[1])
          res$upper[idx, k] <- apply(S_draws, 1, quantile, probs[2])
        }
      }
    }
    colnames(res$mean) <- as.character(times)
    if (quantiles)
      colnames(res$lower) <- colnames(res$upper) <- as.character(times)
  }

  if (!is.null(rownames(X_new))) {
    setnm <- if (type == "survival") `rownames<-` else `names<-`
    for (nm in intersect(c("mean", "lower", "upper"), names(res)))
      res[[nm]] <- setnm(res[[nm]], rownames(X_new))
  }

  return(res)
}
