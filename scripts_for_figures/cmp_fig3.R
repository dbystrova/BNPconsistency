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
data <- loadRData("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_20000.RData")

## number of mixture components
K_nc <- 10
## number of iterations, M without burnin
alpha_1 = 0.01
alpha_2 = 0.5
alpha_3 = 0.9
M_it <- 10000
burnin_ <- 2000

pk_n_001 = MCMC_function(data, e0=alpha_1, K=K_nc, M=M_it, burnin=burnin_) 
pk_n_01 = MCMC_function(data, e0=alpha_2, K=K_nc, M=M_it, burnin=burnin_) 
pk_n_09 = MCMC_function(data, e0=alpha_3, K=K_nc, M=M_it, burnin=burnin_) 


df_= data.frame(K= 1:K_nc, 
                Pkn_1 = pk_n_001$p_k,
                Pkn_2 = pk_n_01$p_k,
                Pkn_3 = pk_n_09$p_k)%>% gather(Process_type, density, Pkn_1:Pkn_3)

df_$alpha = c(rep(alpha_1,K_nc),rep(alpha_2,K_nc),rep(alpha_3,K_nc))  
df_$N = rep(dim(data$y)[1],length(df_$K))
df_$Rh =  c(rep(pk_n_001$ll_rhat,length( pk_n_001$p_k)), rep(pk_n_01$ll_rhat,length( pk_n_001$p_k)),rep(pk_n_09$ll_rhat,length( pk_n_001$p_k)))

save(df_, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig3.RData")
