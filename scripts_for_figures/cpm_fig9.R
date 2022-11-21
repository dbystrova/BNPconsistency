#rm(list=ls()) 
#setwd("~/Documents/GitHub/BNPconsistency/scripts_for_figures")
## read sources
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Random_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Estimation_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Identification_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Gibbs_sampling_function.R")

require(tidyr)
require(e1071)
require(MASS)
require(MCMCpack)
require(mvtnorm)
require(Runuran)
require(flexclust)
library(cowplot)
library(ggplot2)
#---------- B) Specification of the simulation and prior parameters -----------------------------------------------

ds_list<- c("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_20.RData","~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_200.RData",
            "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_2000.RData","~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_20000.RData")



alpha_list<- c(3.2,1.24, 0.81, 0.6)

comparison_n_2<- function(ds_list,K_, M_it, nburn, alpha_l){
  pk<- list()
  N<- c()
  R_h <- c()
  for (i in 1:length(ds_list)){
    data =  loadRData(ds_list[i])
    pk[[i]] <- MCMC_function(data, e0=alpha_l[i]/K_, K=K_, M=M_it, burnin=nburn)
    N[i]<- dim(data$y)[1]
  }
  df_ = tibble(K= 1:K_)
  df2_ = tibble(K= 1:K_)
  df3_ = tibble(K= 1:K_)
  df4_ = tibble(K= 1:K_)
  for (j in 1:length(ds_list)){
    name_ <- paste("Pkn_", j, sep = "")
    name2_ <- paste("Rh_", j, sep = "")
    name3_ <- paste("N_", j, sep = "")
    name4_ <- paste("Alpha_", j, sep = "")
    df_[,name_]<- pk[[j]]$p_k
    df2_[,name2_]<- rep(pk[[j]]$ll_rhat,length(pk[[j]]$p_k))
    df3_[,name3_]<- rep(N[j],length(pk[[j]]$p_k))
    df4_[,name4_]<- rep(alpha_l[j]/K_,length(pk[[j]]$p_k))
  }
  df = df_%>% gather(Process_type, density,  paste("Pkn_", 1, sep = ""):paste("Pkn_", length(ds_list), sep = ""))
  df2 = df2_%>% gather(Rh, Rh_val,  paste("Rh_", 1, sep = ""):paste("Rh_", length(ds_list), sep = ""))
  df3 = df3_%>% gather(N_, N_val,  paste("N_", 1, sep = ""):paste("N_", length(ds_list), sep = ""))
  df4 = df4_%>% gather(Alpha_, Alpha_val,  paste("Alpha_", 1, sep = ""):paste("Alpha_", length(ds_list), sep = ""))
 # df$alpha = c((alpha[1],dim(df_)[1]*length(ds_list))) 
  df$Rh = df2$Rh_val
  df$N =df3$N_val
  df$Al =df4$Alpha_val
  
  df_l_ = tibble(K= 1:((M_it -nburn)*2))
  df2_l_ = tibble(K= 1:((M_it - nburn)*2))
  df3_l_ = tibble(K= 1:((M_it - nburn)*2))
  df4_l_ = tibble(K= 1:((M_it - nburn)*2))
  for (j in 1:length(ds_list)){
    name_ <- paste("P_", j, sep = "")
    name2_ <- paste("Rh_", j, sep = "")
    name3_ <- paste("N_", j, sep = "")
    name4_ <- paste("Alpha_", j, sep = "")
    df_l_[,name_]<- pk[[j]]$p_k_all
    df2_l_[,name2_]<- rep(pk[[j]]$ll_rhat,length(pk[[j]]$p_k_all))
    df3_l_[,name3_]<- rep(N[j],length(pk[[j]]$p_k_all))
    df4_l_[,name4_]<- rep(alpha_list[j],length(pk[[j]]$p_k_all))
  }
  df_l = df_l_%>% gather(Process_type, density,  paste("P_", 1, sep = ""):paste("P_", length(ds_list), sep = ""))
  df_l2 = df2_l_%>% gather(Rh, Rh_val,  paste("Rh_", 1, sep = ""):paste("Rh_", length(ds_list), sep = ""))
  df_l3 = df3_l_%>% gather(N_, N_val,  paste("N_", 1, sep = ""):paste("N_", length(ds_list), sep = ""))
  df_l4 = df4_l_%>% gather(Alpha_, Alpha_val,  paste("Alpha_", 1, sep = ""):paste("Alpha_", length(ds_list), sep = ""))
 # df_l$alpha = c(rep(alpha,(M_it -nburn)*2*length(ds_list))) 
  df_l$Rh = df_l2$Rh_val
  df_l$N =df_l3$N_val
  df_l$Al =df_l4$Alpha_val
  return(list(line = df, hist =df_l))
}




df_9<-comparison_n_2(ds_list,K_ = 10, M_it= 10000 , nburn = 2000,alpha_l = alpha_list)
save(df_9, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig9.RData")
