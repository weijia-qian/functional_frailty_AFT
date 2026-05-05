getPCA <- function(data, nbasis, gp)
{
  
  n <- dim(data)[1]
  p <- dim(data)[2]
  dimnames(data) = list(as.character(1:n), as.character(1:p))
  bs_basis <- create.bspline.basis(rangeval = c(gp[1], gp[p]), nbasis = nbasis)
  inp_mat <- inprod(bs_basis, bs_basis)
  sinp_mat <- expm:::sqrtm(inp_mat)
  evalbase = eval.basis(gp, bs_basis)
  fdobj <- fdPar(bs_basis, int2Lfd(2), lambda=0)
  pcaobj <- smooth.basisPar(gp, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd
  
  mean_coef <- apply(t(pcaobj$coefs), 2, mean)
  sdata <- scale(t(pcaobj$coefs), scale = FALSE)
  new.data <- sdata %*% sinp_mat
  dcov <- cov(new.data)
  d.eigen <- eigen(dcov)
  var_prop <- cumsum(d.eigen$values)/sum(d.eigen$values)
  ncomp <- which(var_prop>0.85)[1]
  loads <- d.eigen$vectors[,1:ncomp]
  PCs <- solve(sinp_mat) %*% loads
  colnames(PCs) = 1:ncomp
  for(i in 1:ncomp)
    colnames(PCs)[i] = paste("PC", i, sep = "")
  PCAcoef <- fd(PCs, bs_basis)
  mean_coef <- fd(as.vector(mean_coef), bs_basis)
  pcaobj2 <- pcaobj
  pcaobj2$coefs <- t(sdata)
  PCAscore <- inprod(pcaobj2, PCAcoef)
  colnames(PCAscore) = 1:ncomp
  for(i in 1:ncomp)
    colnames(PCAscore)[i] = paste("V", i, sep = "")
  
  return(list(PCAcoef = PCAcoef, PCAscore = PCAscore, meanScore = mean_coef,
              bs_basis = bs_basis, evalbase = evalbase, ncomp = ncomp, gp = gp))
}

fpcr <- function(time, event, group = NULL, X, Z, weights = NULL, nb, gp)
{
  
  n <- dim(X)[1]
  j <- dim(X)[2]
  
  fpca <- getPCA(data = X, nbasis = nb, gp = gp)
  fscore <- fpca$PCAscore
  evalbase <- fpca$evalbase
  PCAcoef <- fpca$PCAcoef
  
  model_mat <- cbind(fscore, Z)
  if(is.null(group)){
    id <- 1:n
  }else{
    id <- group
  }
  
  # FIX: Create data.frame directly without cbind()
  model_matrix <- data.frame(time = time, event = event, id = id, model_mat)
  
  for(i in 4:dim(model_matrix)[2]) {
    colnames(model_matrix)[i] = paste("V", (i-3), sep = "")
  }
  
  var_name <- paste(c(colnames(model_matrix)[-(1:3)], "(1|id)"), collapse = "+")
  
  cox_formula <- as.formula(paste("Surv(time, event)~", var_name, sep = ""))
  if(is.null(weights)){
    cox_model <- coxme(cox_formula, data = model_matrix)
  }else{
    cox_model <- coxme(cox_formula, data = model_matrix, weights = weights,
                       control = coxme.control(
                         eps = 1e-05, toler.chol = .Machine$double.eps^0.75,
                         iter.max = 5000, sparse.calc = NULL,
                         optpar = list(method = "CG",
                                       control = list(reltol = 1e-3, maxit = 5000))))
  }
  
  cox_result <- summary(cox_model)
  cox_coefs <- cox_result$coefficients
  coef_f <- cox_coefs[,1]
  
  conc <- 1 - concordance(Surv(time, event) ~
                            cox_model$linear.predictor, model_matrix)$concordance
  
  bhat_t <- evalbase %*% (PCAcoef$coefs %*% as.matrix(coef_f[1:dim(fscore)[2]]))
  gamma_hat <- coef_f[-(1:dim(fscore)[2])]
  
  return(list(bhat = bhat_t, gammahat = gamma_hat, concordance = conc,
              model = cox_model, fpca = fpca))
  
}

getPCA_test <- function(object, data)
{
  
  bs_basis <- object$bs_basis
  PCAcoef <- object$PCAcoef
  gp <- object$gp
  mean.tr <- c(object$meanScore$coefs)
  n <- dim(data)[1]
  p <- dim(data)[2]
  dimnames(data) = list(as.character(1:n), as.character(1:p))
  pcaobj <- smooth.basisPar(gp, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd
  
  sdata <- scale(t(pcaobj$coefs), center = mean.tr, scale = FALSE)
  pcaobj2 <- pcaobj
  pcaobj2$coefs <- t(sdata)
  PCAscore_test = inprod(pcaobj2, PCAcoef)
  colnames(PCAscore_test) = 1:dim(PCAcoef$coefs)[2]
  for(i in 1:dim(PCAcoef$coefs)[2])
    colnames(PCAscore_test)[i] = paste("Score", i, sep = "")
  
  return(PCAscore.test = PCAscore_test)
}

#' Calculate Survival Probabilities from Cox Linear Predictors
#' 
#' @param eta_train Vector of linear predictors from the training set
#' @param time_train Vector of survival times from the training set
#' @param event_train Vector of event indicators from the training set (1 = event, 0 = censored)
#' @param eta_test Vector of linear predictors for the test set
#' @param tgrid Vector of time points to predict survival probabilities
get_breslow_survival <- function(eta_train, time_train, event_train, eta_test, tgrid) {
  
  # 1. Order training data by time
  ord <- order(time_train)
  t_train  <- time_train[ord]
  e_train  <- event_train[ord]
  lp_train <- eta_train[ord]
  
  # Risk scores: exp(eta)
  risk_train <- exp(lp_train)
  risk_test  <- exp(eta_test)
  
  # 2. Identify unique event times
  unique_events <- unique(t_train[e_train == 1])
  
  # 3. Calculate Breslow baseline hazard increments
  h0 <- numeric(length(unique_events))
  for (i in seq_along(unique_events)) {
    t_ev <- unique_events[i]
    
    # Number of events exactly at t_ev
    d_i <- sum(t_train == t_ev & e_train == 1)
    
    # Sum of risk scores for all individuals still at risk (T_i >= t_ev)
    sum_risk <- sum(risk_train[t_train >= t_ev])
    
    h0[i] <- d_i / sum_risk
  }
  
  # 4. Cumulative baseline hazard: H0(t)
  H0 <- cumsum(h0)
  
  # 5. Evaluate H0 at the requested tgrid using a step function
  H0_tgrid <- numeric(length(tgrid))
  for (k in seq_along(tgrid)) {
    t_star <- tgrid[k]
    valid_idx <- which(unique_events <= t_star)
    
    if (length(valid_idx) == 0) {
      H0_tgrid[k] <- 0  # No hazard accumulated yet
    } else {
      H0_tgrid[k] <- H0[max(valid_idx)]
    }
  }
  
  # 6. Compute S(t|X) for the test set (N_test x length(tgrid) matrix)
  S_mat <- matrix(NA, nrow = length(eta_test), ncol = length(tgrid))
  for (k in seq_along(tgrid)) {
    S_mat[, k] <- exp(-H0_tgrid[k] * risk_test)
  }
  
  return(S_mat)
}

fpcr_predict <- function(object, time, event, group = NULL, X, Z, nb, gp)
{
  
  n <- dim(X)[1]
  j <- dim(X)[2]
  
  fpca <- object$fpca
  fscore <- getPCA_test(fpca, X)
  evalbase <- fpca$evalbase
  PCAcoef <- fpca$PCAcoef
  
  model_mat <- cbind(fscore, Z)
  if(is.null(group)){
    id <- 1:n
  }else{
    id <- group
  }
  
  # MINIMAL FIX 1: Remove cbind() here so model_mat stays strictly numeric
  model_matrix <- data.frame(time=time, event=event, id=id, model_mat)
  
  for(i in 4:dim(model_matrix)[2])
    colnames(model_matrix)[i] = paste("V", (i-3), sep = "")
  
  colnames(model_matrix)[1] <- "time"
  
  var_name <- paste(c(colnames(model_matrix)[-(1:3)], "(1|id)"), collapse = "+")
  
  cox_model <- object$model
  
  # MINIMAL FIX 2: Replace predict_coxme with inline matrix multiplication
  preds <- as.numeric(as.matrix(model_matrix[, -(1:3)]) %*% coxme::fixef(cox_model))
  
  conc <- 1 - concordance(Surv(time, event) ~ preds, model_matrix)$concordance
  
  return(list(predictions = preds, concordance = conc))
}
