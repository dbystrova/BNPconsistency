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

####### generate data number of generating components
#number of components
K_g <- 3
# mixture weights
eta <- c(0.4,0.3,0.3)
mu <- cbind(c(0.8,0.8), c(0.8,-0.8), c(-0.8,0.8))
## covariance matrices:
r <- length(mu[, 1])
sigma <- array(0, dim = c(r, r, K_g))
for (k in 1:K_g) {
  sigma[, , k] <- diag(c(rep(0.05, r)))
}

####### generate data number of generating components

## number of observations
#N <- 1000
## create dataset
N=500
date_500 <- MultVar_RandomNormal(N=500, eta, mu, sigma, z = rep(FALSE, N))
save(date_500, file = "sim_data/GM_3_500.RData")

N=1500
date_1500 <- MultVar_RandomNormal(N=1500, eta, mu, sigma, z = rep(FALSE, N))
save(date_1500, file = "sim_data/GM_3_1500.RData")

N=5000

date_5000 <- MultVar_RandomNormal(N=5000, eta, mu, sigma, z = rep(FALSE, N))
save(date_5000, file = "sim_data/GM_3_5000.RData")

N=10000
date_10000 <- MultVar_RandomNormal(N=10000, eta, mu, sigma, z = rep(FALSE, N))
save(date_10000, file = "sim_data/GM_3_10000.RData")
