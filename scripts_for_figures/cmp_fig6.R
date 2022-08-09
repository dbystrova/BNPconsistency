#rm(list=ls()) 
#setwd("~/Documents/GitHub/BNPconsistency/scripts_for_figures")
## read sources
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Random_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Estimation_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Identification_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Gibbs_sampling_function.R")

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
data_500 <- loadRData("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_500.RData")
data_1500 <- loadRData("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_1500.RData")
data_5000 <- loadRData("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_5000.RData")

## number of mixture components
K_nc <- 10
## number of iterations, M without burnin
alpha_1 = 0.05
#alpha_2 = 0.1
#alpha_3 = 0.9
M_it <- 10000
burnin_ <- 2000

pk_n_500 = MCMC_function(data_500, e0=alpha_1, K=K_nc, M=M_it, burnin=burnin_) 
pk_n_1500 = MCMC_function(data_1500, e0=alpha_1, K=K_nc, M=M_it, burnin=burnin_) 
pk_n_5000 = MCMC_function(data_5000, e0=alpha_1, K=K_nc, M=M_it, burnin=burnin_) 


df_= data.frame(K= 1:K_nc, 
                Pkn_1 = pk_n_500,
                Pkn_2 = pk_n_1500,
                Pkn_3 = pk_n_5000)%>% gather(Process_type, density, Pkn_1:Pkn_3)

df_$alpha = c(rep(alpha_1,K_nc),rep(alpha_1,K_nc),rep(alpha_1,K_nc))  
df_$N = c(rep(dim(data_500$y)[1],K_nc),rep(dim(data_1500$y)[1],K_nc),rep(dim(data_5000$y)[1],K_nc)) 
save(df_, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig6.RData")
