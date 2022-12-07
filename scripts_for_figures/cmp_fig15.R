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

comparison_n_3<- function(ds_list,K_, M_it, nburn){
  pk<- list()
  N<- c()
  R_h <- c()
  W_non_sorted <- list()
  W <- list()
  Mu_mat <- list()
  S_mat<- list()
  for (i in 1:length(ds_list)){
    data =  loadRData(ds_list[i])
    pk[[i]] <- MCMC_function(data, K=K_, M=M_it, burnin=nburn, priorOnE0 = TRUE)
    N[i]<- dim(data$y)[1]
    #Eta_ <- apply(pk[[1]]$Eta , 2, function(x) c(sort(x, decreasing = TRUE)))
    Eta_<- matrix(NA, nrow =dim(pk[[i]]$Eta)[1],ncol =  dim(pk[[i]]$Eta)[2] )
    for (j in 1:dim(pk[[i]]$Eta)[1]){ Eta_[j,] <- sort(pk[[i]]$Eta[j,],decreasing = TRUE)}
    W[[i]] <- Eta_
    W_non_sorted[[i]] <- pk[[i]]$Eta
    Mu_mat[[i]] <- pk[[i]]$Mu
    S_mat[[i]] <-pk[[i]]$Sigma
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
  df4_post <- cbind(df4,t(W_df))
  
  df4_post_ <- gather(df4_post, key = "it",value="weights", 4: dim(df4_post)[2])
  
  #df$alpha = c(rep(alpha,dim(df_)[1]*length(ds_list))) 
  df$Rh = df2$Rh_val
  df$N =df3$N_val
  
  df_l_ = tibble(K= 1:((M_it -nburn)*2))
  df2_l_ = tibble(K= 1:((M_it - nburn)*2))
  df3_l_ = tibble(K= 1:((M_it - nburn)*2))
  for (j in 1:length(ds_list)){
    name_ <- paste("P_", j, sep = "")
    name2_ <- paste("Rh_", j, sep = "")
    name3_ <- paste("N_", j, sep = "")
    df_l_[,name_]<- pk[[j]]$p_k_all
    df2_l_[,name2_]<- rep(pk[[j]]$ll_rhat,length(pk[[j]]$p_k_all))
    df3_l_[,name3_]<- rep(N[j],length(pk[[j]]$p_k_all))
  }
  df_l = df_l_%>% gather(Process_type, density,  paste("P_", 1, sep = ""):paste("P_", length(ds_list), sep = ""))
  df_l2 = df2_l_%>% gather(Rh, Rh_val,  paste("Rh_", 1, sep = ""):paste("Rh_", length(ds_list), sep = ""))
  df_l3 = df3_l_%>% gather(N_, N_val,  paste("N_", 1, sep = ""):paste("N_", length(ds_list), sep = ""))
  #df_l$alpha = c(rep(alpha,(M_it -nburn)*2*length(ds_list))) 
  df_l$Rh = df_l2$Rh_val
  df_l$N =df_l3$N_val
  return(list(line = df, hist =df_l, weights = df4_post_, eta = W_non_sorted, mu = Mu_mat, sig = S_mat))
}

M_it <- 100
burnin_ <- 20

ds_list<- c("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_20.RData","~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_200.RData",
            "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_2000.RData","~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_20000.RData")



df_15<-comparison_n_3(ds_list,K_ = 10, M_it= 10000 , nburn = 2000)
save(df_15, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig15.RData")
