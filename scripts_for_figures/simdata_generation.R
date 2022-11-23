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
eta <- c(0.4,0.3,0.2)
mu <- cbind(c(0.8,0.8), c(0.8,-0.8), c(-0.8,0.8))
## covariance matrices:
r <- length(mu[, 1])
sigma <- array(0, dim = c(r, r, K_g))
for (k in 1:K_g) {
  sigma[, , k] <- diag(c(rep(0.05, r)))
}

####### generate data number of generating components

N=20000
date_20000 <- MultVar_RandomNormal(N=20000, eta, mu, sigma, z = rep(FALSE, N),seed=12345)
save(date_20000, file = "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_20000.RData")
write.csv(date_20000$y,"~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_20000.csv", row.names = FALSE)

date_5000 =list()
date_5000$eta = date_20000$eta 
date_5000$mu = date_20000$mu 
date_5000$sigma = date_20000$sigma 

set.seed(1234)
date_5000$y <- date_20000$y[sample(1:nrow(date_20000$y), 5000,replace=FALSE),]
date_5000$z <-rep(FALSE, 5000)
save(date_5000, file = "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_5000.RData")
write.csv(date_5000$y,"~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_5000.csv", row.names = TRUE)

date_2000 =list()
date_2000$eta = date_5000$eta 
date_2000$mu = date_5000$mu 
date_2000$sigma = date_5000$sigma 

date_2000$y <- date_5000$y[sample(1:nrow(date_5000$y), 2000,replace=FALSE),]
date_2000$z <-rep(FALSE, 2000)
save(date_2000, file = "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_2000.RData")
write.csv(date_2000$y,"~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_2000.csv", row.names = TRUE)

date_200 =list()
date_200$eta = date_2000$eta 
date_200$mu = date_2000$mu 
date_200$sigma = date_2000$sigma 

date_200$y <- date_2000$y[sample(1:nrow(date_2000$y), 200,replace=FALSE),]
date_200$z <-rep(FALSE, 200)
save(date_200, file = "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_200.RData")
write.csv(date_200$y,"~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_200.csv", row.names = TRUE)


date_20 =list()
date_20$eta = date_200$eta 
date_20$mu = date_200$mu 
date_20$sigma = date_200$sigma 

date_20$y <- date_200$y[sample(1:nrow(date_200$y), 20,replace=FALSE),]
date_20$z <-rep(FALSE, 20)
save(date_20, file = "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_20.RData")
write.csv(date_20$y,"~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_20.csv", row.names = FALSE)

#overlap <- subset(date_1500$y,date_1500$y[,1] %in% intersect(date_1500$y[,1],date_10000$y[,1]))
