bs_coef = c(0, -0.6, -1.2, 0.6, -0.5, 1, 0.5, 0)
bs_coef <- c(0.0, -0.4, 0.6, 1.4, 1.6, 1.5, 1.3, 1.0, 0.4, -0.2, -0.5, 0.0)
# Option A shape, scaled to have similar effect size to your "simple" + "complex"
bs_coef <- c(0.0, -0.6, -0.4, 0.2, 0.7, 0.6, 0.1, -0.3, -0.5, -0.1)
bs_coef <- c(
  0.0,   0.9,  1.1,  0.4, -0.6,
  -1.2,  -0.8,  0.2,  1.0,  1.3,
  0.7,  -0.2, -0.7, -0.3,  0.0
)

bs_coef = c(0.0, 1.7, -2.1, 2.5, -2.3, 2.1, -1.9, 1.5, -1.25, 1.1, -0.85, 0.6, -0.4, 0.2, 0.0)

k <- length(bs_coef)
B <- splines::bs(s_grid, df = k, intercept = TRUE) # ns x k
beta <- B %*% bs_coef
plot(s_grid, beta, type="l", main=paste0("beta(s), k=", k), xlab="s", ylab="beta(s)")
abline(h=0, lty=2)


beta_true <- sim_data$coefficients$beta
beta_true_int <- mean(beta_true)
beta_est <- fit$beta_mean
beta_est_int <- mean(beta_est)
beta_true_int; beta_est_int

matplot(t(sim_curves[1:15,]), type = "l", xlab = "Function index", ylab = "X(s)")
matplot(t(sim_data$data$X[1:10,]), type = "l", xlab = "Function index", ylab = "X(s)")

simulate_X <- function(N, s_grid, K0 = 20, alpha = 0.5) {
  
  nS <- length(s_grid)
  
  # ---- Construct Fourier basis ----
  Phi <- matrix(0, nS, K0)
  
  Phi[,1] <- 1
  k_idx <- 1
  
  for (k in 1:floor((K0-1)/2)) {
    Phi[,2*k]   <- sqrt(2) * cos(2*pi*k*s_grid)
    Phi[,2*k+1] <- sqrt(2) * sin(2*pi*k*s_grid)
  }
  
  # ---- Eigenvalue decay ----
  lambda <- (1:K0)^(-alpha)
  
  # ---- Random scores ----
  Xi <- matrix(rnorm(N*K0), N, K0)
  
  # ---- Generate curves ----
  X <- Xi %*% diag(sqrt(lambda)) %*% t(Phi)
  
  return(X)
}

X <- simulate_X(N=200, s_grid, K0=15, alpha=0)

sv <- svd(X)$d
p <- sv / sum(sv)
rank_eff <- exp(-sum(p * log(p)))

rank_eff