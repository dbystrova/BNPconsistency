
## read sources
source("Code_SP_Mix/Random_SpMix.R")
source("Code_SP_Mix/Estimation_SpMix.R")
source("Code_SP_Mix/Identification_SpMix.R")

require(e1071)
require(mclust)
require(MASS)
require(bayesm)
require(MCMCpack)
require(mvtnorm)
require(Runuran)
require(flexclust)
library(gridExtra)
library(cowplot)
library(ggplot2)
#---------- B) Specification of the simulation and prior parameters -----------------------------------------------

MCMC_function <- function(data, e0=0.01, K, M, burnin) {
  y <- as.matrix(data$y)
  Mmax <- M + burnin
  
  ## read dimensions of data:
  r <- length(y[1, ])
  N <- length(y[, 1])
  
  ## Dirichlet parameter for the mixture weights
  ## variance of the normal proposal for the MH step for estimating e0
  c_proposal <- 0.8
  
  R <- apply(y, 2, function(x) diff(range(x)))
  
  ## prior on Sigma_k
  c0 <- 2.5 + (r - 1)/2
  C0 <- 0.75 * cov(y)
  g0 <- 0.5 + (r - 1)/2
  G0 <- 100 * g0/c0 * diag((1/R^2))
  
  ## prior on mu
  b0 <- apply(y, 2, median)
  B_0 <- rep(1, r)  #initial values for lambda are 1
  B0 <- diag((R^2) * B_0)
  nu <- 0.5
  
  ## initial values for parameters to be estimated:
  eta_0 <- rep(1/K, K)
  sigma_0 <- array(0, dim = c(r, r, K))
  for (k in 1:K) {
    sigma_0[, , k] <- C0
  }
  C0_0 <- C0
  
  ## initial classification
  groups <- K
  cl_y <- kmeans(y, centers = groups, nstart = 30)
  S_0 <- cl_y$cluster
  mu_0 <- cbind(t(cl_y$centers))
  
  
  ## generate matrices for saving the results:
  Eta <- matrix(0, M, K)
  Mu <- array(0, dim = c(M, r, K))
  B <- matrix(0, M, r)
  Eta_Matrix_FS <- matrix(0, r, K)
  Sigma_Matrix_FS <- array(0, dim = c(r, r, K))
  Mu_Matrix_FS <- array(0, dim = c(r, K))
  
  
  
  #---------- C) Gibbs sampling from the posterior -----------------------------------------------
  print(B_0)
  ################ call MCMC procedure
  estGibbs <- MultVar_NormMixt_Gibbs_IndPriorNormalgamma(y, S_0, mu_0, sigma_0, eta_0, e0, c0, C0_0, 
                                                         g0, G0, b0, B0, nu, B_0, M, burnin, c_proposal, priorOnE0 = FALSE, lambda = FALSE)
  
  Mu <- estGibbs$Mu
  Sigma <- estGibbs$Sigma
  Eta <- estGibbs$Eta
  S_alt_matrix <- estGibbs$S_alt_matrix
  B <- estGibbs$B
  e0_vector <- estGibbs$e0_vector
  acc_rate <- estGibbs$acc_rate
  Nk_matrix_alt <- estGibbs$Nk_matrix_alt
  nonnormpost_mode_list <- estGibbs$nonnormpost_mode_list
  
  ##### number of nonempty components (nne_gr), after burnin:
  nne_gr <- apply(Nk_matrix_alt, 1, function(x) sum(x != 0))
  table(nne_gr)
  #plot(1:length(nne_gr), nne_gr, type = "l", ylim = c(1, max(nne_gr)), 
  #     xlab = paste("iteration"), main = "number of non-empty components")
  
  #---------- E) Identification of the mixture model -----------------------------------------------
  ##### estimating the number of nonempty components
  K0_vector <- rowSums(Nk_matrix_alt != 0)  #vector with number of non-empty groups of each iteration
  p_K0 <- tabulate(K0_vector, K)
  p_K0
  #par(mfrow = c(1, 1))
  #barplot(p_K0, names = 1:K, xlab = "number of non-empty groups K0", col = "green", ylab = "freq")
  K0 <- which.max(p_K0)
  K0  #mode K0 is the estimator for K_true
  M0 <- sum(K0_vector == K0)
  M0  #M0 draws have exactly K0 non-empty groups
return(p_K0)
}
