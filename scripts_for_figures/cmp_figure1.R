#rm(list=ls()) 
#setwd("~/Documents/Code Consistency/Code_SpMix")
## read sources
source("Random_SpMix.R")
source("Estimation_SpMix.R")
source("Identification_SpMix.R")
source("Gibbs_sampling_function.R")

require(tidyr)
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


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
data <- loadRData("sim_data/GM_3_500.RData")

## number of mixture components
K_nc <- 10
## number of iterations, M without burnin
alpha_1 = 0.01
alpha_2 = 0.1
alpha_3 = 0.9
M_it <- 5000
burnin_ <- 2000

pk_n_001 = MCMC_function(data, e0=alpha_1, K=K_nc, M=M_it, burnin=burnin_) 
pk_n_01 = MCMC_function(data, e0=alpha_2, K=K_nc, M=M_it, burnin=burnin_) 
pk_n_09 = MCMC_function(data, e0=alpha_3, K=K_nc, M=M_it, burnin=burnin_) 


df_= data.frame(K= 1:K, 
                Pkn_1 = pk_n_001,
                Pkn_2 = pk_n_01,
                Pkn_3 = pk_n_09)%>% gather(Process_type, density, Pkn_1:Pkn_3)

df_$alpha = c(rep(alpha_1,K),rep(alpha_2,K),rep(alpha_3,K))                
save(df_, file = "sim_data/cmp_fig1.RData")
