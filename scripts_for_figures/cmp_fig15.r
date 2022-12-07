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

ds_list<- c("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_20.RData","~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_200.RData",
            "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_2000.RData","~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_20000.RData")

comparison_n_3<- function(ds_list,K_, M_it, nburn){
  pk<- list()
  N<- c()
  R_h <- c()
  W_non_sorted <- list()
  W <- list()
  Mu_mat <- list()
  S_mat<- list()
  E0_list<- list()
  for (i in 1:length(ds_list)){
    data =  loadRData(ds_list[i])
    pk[[i]] <- MCMC_function_ra(data, K=K_, M=M_it, burnin=nburn)
    N[i]<- dim(data$y)[1]
    Eta_<- matrix(NA, nrow =dim(pk[[i]]$Eta)[1],ncol =  dim(pk[[i]]$Eta)[2] )
    for (j in 1:dim(pk[[i]]$Eta)[1]){ Eta_[j,] <- sort(pk[[i]]$Eta[j,],decreasing = TRUE)}
    W[[i]] <- Eta_
    W_non_sorted[[i]] <- pk[[i]]$Eta
    Mu_mat[[i]] <- pk[[i]]$Mu
    S_mat[[i]] <-pk[[i]]$Sigma
    E0_list[[i]] <-pk[[i]]$e0
  }
  df_ = tibble(K= 1:K_)
  df2_ = tibble(K= 1:K_)
  df3_ = tibble(K= 1:K_)
  df4_ = tibble(K= 1:K_)
  for (j in 1:length(ds_list)){
    name_ <- paste("Pkn_", j, sep = "")
    name2_ <- paste("Rh_", j, sep = "")
    name3_ <- paste("N_", j, sep = "")
    name4_ <- paste("W_", j, sep = "")
    df_[,name_]<- pk[[j]]$p_k
    df2_[,name2_]<- rep(pk[[j]]$ll_rhat,length(pk[[j]]$p_k))
    df3_[,name3_]<- rep(N[j],length(pk[[j]]$p_k))
    df4_[,name4_]<- rep(N[j],length(pk[[j]]$p_k))
  }
  df = df_%>% gather(Process_type, density,  paste("Pkn_", 1, sep = ""):paste("Pkn_", length(ds_list), sep = ""))
  df2 = df2_%>% gather(Rh, Rh_val,  paste("Rh_", 1, sep = ""):paste("Rh_", length(ds_list), sep = ""))
  df3 = df3_%>% gather(N_, N_val,  paste("N_", 1, sep = ""):paste("N_", length(ds_list), sep = ""))
  df4 = df4_%>% gather(W_, W_val,  paste("W_", 1, sep = ""):paste("W_", length(ds_list), sep = ""))
 
  W_df <- do.call(cbind, W)
  df5_post <- cbind(df5,t(W_df))
  
  df5_post_ <- gather(df5_post, key = "it",value="weights", 4: dim(df5_post)[2])

  
  df$Rh = df2$Rh_val
  df$N =df3$N_val

  
  df_l_ = tibble(it= 1:((M_it -nburn)*2))
  df2_l_ = tibble(it= 1:((M_it - nburn)*2))
  df3_l_ = tibble(it= 1:((M_it - nburn)*2))
  df5_l_ = tibble(it= 1:((M_it - nburn)*2))
  for (j in 1:length(ds_list)){
    name_ <- paste("P_", j, sep = "")
    name2_ <- paste("Rh_", j, sep = "")
    name3_ <- paste("N_", j, sep = "")
    name5_ <- paste("E0_", j, sep = "")
    df_l_[,name_]<- pk[[j]]$p_k_all
    df2_l_[,name2_]<- rep(pk[[j]]$ll_rhat,length(pk[[j]]$p_k_all))
    df3_l_[,name3_]<- rep(N[j],length(pk[[j]]$p_k_all))
    df5_l_[,name5_]<-  E0_list[[i]]
  }
  df_l = df_l_%>% gather(Process_type, density,  paste("P_", 1, sep = ""):paste("P_", length(ds_list), sep = ""))
  df_l2 = df2_l_%>% gather(Rh, Rh_val,  paste("Rh_", 1, sep = ""):paste("Rh_", length(ds_list), sep = ""))
  df_l3 = df3_l_%>% gather(N_, N_val,  paste("N_", 1, sep = ""):paste("N_", length(ds_list), sep = ""))
  df_l5 = df5_l_%>% gather(E0_, E0_val,  paste("E0_", 1, sep = ""):paste("E0_", length(ds_list), sep = ""))
  df_l$Rh = df_l2$Rh_val
  df_l$N =df_l3$N_val
  df_l5$N = df_l3$N_val
  return(list(line = df, hist =df_l, weights = df5_post_, mu = Mu_mat, sig = S_mat, e0 = df_l5))
}



df_ = comparison_n_3(ds_list,K_ = 10, M_it = 10000, nburn= 2000)

save(df_, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig15.RData")



