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
eta <- c(0.5,0.3,0.2)
mu <- cbind(c(0.8,0.8), c(0.8,-0.8), c(-0.8,0.8))
## covariance matrices:
r <- length(mu[, 1])
sigma <- array(0, dim = c(r, r, K_g))
for (k in 1:K_g) {
  sigma[, , k] <- diag(c(rep(0.05, r)))
}

gen_subsample<- function(data, N, seed_ = 1234){
  set.seed(seed_)
  data_sub =list()
  data_sub$eta = data$eta 
  data_sub$mu = data$mu 
  data_sub$sigma = data$sigma 
  ind_ = sample(1:nrow(data$y), N,replace=FALSE)
  data_sub$y <-  data$y[ind_,]
  data_sub$z <-data$z[ind_]
  save(data_sub, file =paste0("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_",N,".RData"))
  write.csv(data_sub$y,paste0("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_",N,".csv"), row.names = FALSE)
  return(data_sub)
}


####### generate data number of generating components
N=40000
data_40000 <- MultVar_RandomNormal(N=40000, eta, mu, sigma, z = rep(FALSE, N),seed=12345)
save(data_40000, file = "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_40000.RData")
write.csv(data_40000$y,"~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_40000.csv", row.names = FALSE)

data_20000 = gen_subsample(data= data_40000, N= 20000, seed_ = 1234)

data_5000 = gen_subsample(data= data_20000, N= 5000, seed_ = 1234)

data_2000= gen_subsample(data= data_5000, N= 2000, seed_ = 1234)

data_200= gen_subsample(data= data_2000, N= 200, seed_ = 1234)

data_20= gen_subsample(data= data_200, N= 20, seed_ = 1234)

