## read sources
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Random_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Estimation_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Identification_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Gibbs_sampling_function.R")

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
#N=500
#date_500 <- MultVar_RandomNormal(N=500, eta, mu, sigma, z = rep(FALSE, N))
#save(date_500, file = "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_500.RData")

#N=1500
#date_1500 <- MultVar_RandomNormal(N=1500, eta, mu, sigma, z = rep(FALSE, N))
#save(date_1500, file = "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_1500.RData")

#N=5000
#date_5000 <- MultVar_RandomNormal(N=5000, eta, mu, sigma, z = rep(FALSE, N))
#save(date_5000, file = "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_5000.RData")

N=10000
date_10000 <- MultVar_RandomNormal(N=10000, eta, mu, sigma, z = rep(FALSE, N))
save(date_10000, file = "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_10000.RData")

date_5000 =date_10000 
date_5000$y <- date_10000$y[sample(1:nrow(date_10000$y), 5000,replace=FALSE),]
date_5000$z <-rep(FALSE, 5000)
save(date_5000, file = "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_5000.RData")


date_1500 =date_5000 
date_1500$y <- date_5000$y[sample(1:nrow(date_5000$y), 1500,replace=FALSE),]
date_1500$z <-rep(FALSE, 1500)
save(date_1500, file = "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_1500.RData")


date_500 =date_1500 
date_500$y <- date_1500$y[sample(1:nrow(date_1500$y), 500,replace=FALSE),]
date_500$z <-rep(FALSE, 500)
save(date_500, file = "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_500.RData")


#overlap <- subset(date_1500$y,date_1500$y[,1] %in% intersect(date_1500$y[,1],date_10000$y[,1]))
