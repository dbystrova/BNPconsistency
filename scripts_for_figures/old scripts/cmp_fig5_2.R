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
data_10 <- loadRData("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_10.RData")
data_100 <- loadRData("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_100.RData")
data_1000 <- loadRData("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_1000.RData")
data_10000 <- loadRData("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_10000.RData")

## number of mixture components
K_nc <- 10
## number of iterations, M without burnin
#alpha_1 = 0.01
alpha_2 = 0.5
#alpha_3 = 0.9
#alpha_4 =1
#alpha_5 = 2
M_it <- 5000
burnin_ <- 2000
K_nc <- 9
pk_n_10 = MCMC_function(data_10, e0=alpha_2, K=K_nc, M=M_it, burnin=burnin_) 
K_nc <- 10
pk_n_100 = MCMC_function(data_100, e0=alpha_2, K=K_nc, M=M_it, burnin=burnin_) 
pk_n_1000 = MCMC_function(data_1000, e0=alpha_2, K=K_nc, M=M_it, burnin=burnin_) 
pk_n_10000 = MCMC_function(data_10000, e0=alpha_2, K=K_nc, M=M_it, burnin=burnin_) 


df_= data.frame(K= 1:K_nc, 
                Pkn_1 = c(pk_n_10$p_k,0),
                Pkn_2 = pk_n_100$p_k,
                Pkn_3 = pk_n_1000$p_k,
                Pkn_4 = pk_n_10000$p_k)%>% gather(Process_type, density, Pkn_1:Pkn_4)
df_$alpha = c(rep(alpha_2,K_nc),rep(alpha_2,K_nc),rep(alpha_2,K_nc),rep(alpha_2,K_nc))  
df_$N = c(rep(dim(data_10$y)[1],K_nc),rep(dim(data_100$y)[1],K_nc),rep(dim(data_1000$y)[1],K_nc),rep(dim(data_10000$y)[1],K_nc)) 
df_$Rh =  c(rep(pk_n_10$ll_rhat,K_nc), rep(pk_n_100$ll_rhat,K_nc),rep(pk_n_1000$ll_rhat,K_nc),rep(pk_n_10000$ll_rhat,K_nc))

save(df_, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig5.RData")

df2_= data.frame(Pkn_1 = pk_n_10$p_k_all,
                 Pkn_2 = pk_n_100$p_k_all,
                 Pkn_3 = pk_n_1000$p_k_all,
                 Pkn_4 = pk_n_10000$p_k_all)%>% gather(Process_type, density, Pkn_1:Pkn_4)

N_gr = length( pk_n_10$p_k_all)

df2_$Rh_mu <- c(rep(pk_n_10$m_rh[4],N_gr),rep(pk_n_100$m_rh[4],N_gr),rep(pk_n_1000$m_rh[4],N_gr),rep(pk_n_10000$m_rh[4],N_gr))  
df2_$Rh_s <- c(rep(pk_n_10$s_rh[4],N_gr),rep(pk_n_100$s_rh[4],N_gr),rep(pk_n_1000$s_rh[4],N_gr),rep(pk_n_10000$s_rh[4],N_gr))  
df2_$N = c(rep(dim(data_10$y)[1],N_gr),rep(dim(data_100$y)[1],N_gr),rep(dim(data_1000$y)[1],N_gr),rep(dim(data_10000$y)[1],N_gr)) 
df2_$Rh =  c(rep(pk_n_10$ll_rhat,N_gr), rep(pk_n_100$ll_rhat,N_gr),rep(pk_n_1000$ll_rhat,N_gr),rep(pk_n_10000$ll_rhat,N_gr))
save(df2_, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig5_2.RData")




