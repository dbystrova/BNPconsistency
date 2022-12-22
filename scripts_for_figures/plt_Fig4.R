#rm(list=ls()) 
#setwd("~/Documents/GitHub/BNPconsistency/scripts_for_figures")
## read sources
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
library(latex2exp)
require(tidyr)
library(dplyr)
library(JuliaCall)
library(viridis)
library(bayesm)
library(purrr)
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")

#---------- B) Specification of the simulation and prior parameters -----------------------------------------------

#input_file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig15_2.RData"

simulatuion_function_DPM<- function(nu1=10,nu2= 10*10,ns,Sn,N_tr){
  alpha_s<- rgamma(ns, nu1,nu2)
  alpha_s_mod<- replace(alpha_s, alpha_s< 10^(-10), 10^(-10)) #to avoid small values for alpha, which could lead to inf values in funcDP/PYmo
  p_list<- matrix(NA, nrow =length(alpha_s_mod), ncol=10 )
  julia_library("GibbsTypePriors")
  for (i in 1:length(alpha_s_mod)){
    julia_assign("al", alpha_s_mod[i])
    julia_assign("n", Sn)
    p_list[i,] =  round(julia_eval("Pkn_Dirichlet_mult.(1:10,n, 10, al)"),3)
  }
 return(colMeans(p_list))
}

gelman_diag<- function(c){
  e0_c <- mcmc.list(mcmc(c[1:8000]),mcmc(c[8001:16000]))
  Rhat_al<- gelman.diag(e0_c)$psrf[1]
  return(Rhat_al)
  
}


plt_fig2<-function(input_file, c_v =c(0.1, 0.5, 1, 2), n_list= c(20,200,2000,20000), fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure4/" , a_g = 10, b_g = 10){
  
  # fig_path = "../figures/Figure1/"
  fig_df <- loadRData(input_file)
  fig_df_mut <- fig_df$line%>%group_by(Process_type,N)%>%mutate(pkn =density/sum(density))
  K_= max(fig_df_mut$K)
  julia <- julia_setup()
  julia_library("GibbsTypePriors")
  
  
  
  julia <- julia_setup()
  julia_library("GibbsTypePriors")
   
  df_prior = tibble(K= 1:10, 
                    Pkn_1 = simulatuion_function_DPM(nu1 = a_g, nu2 = K_*b_g,ns = 100, Sn=20, N_tr = 10 ),
                    Pkn_2= simulatuion_function_DPM( nu1 = a_g, nu2 = K_*b_g,ns = 100,  Sn=200, N_tr = 10 ), 
                    Pkn_3 = simulatuion_function_DPM(  nu1 = a_g, nu2 = K_*b_g,ns = 100, Sn=2000, N_tr = 10 ),
                    Pkn_4 = simulatuion_function_DPM(nu1 = a_g, nu2 = K_*b_g,ns = 100, Sn=20000, N_tr = 10 ))%>% gather(Process_type, pkn,Pkn_1:Pkn_4)
  
  df_prior$Type =rep("Prior", dim(df_prior)[1])
  fig_df_mut$Type= rep("Posterior", dim(fig_df_mut)[1])
  
  df_merged = rbind(df_prior,fig_df_mut[, c("K","Process_type", "pkn", "Type")] )
  
  pkn.labs<- c("Pkn_1","Pkn_2","Pkn_3","Pkn_4")
  names(pkn.labs) <- c(paste0("n = ", fig_df_mut$N[1]),paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]),paste0("n = ", fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),paste0("n = ", fig_df_mut$N[(3*max(fig_df_mut$K)+1)]))
  pkn_names <- as_labeller(
    c(`Pkn_1` = paste0("n = ", fig_df_mut$N[1]), `Pkn_2` = paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]),`Pkn_3` = paste0("n = ", fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),`Pkn_4` = paste0("n = ", fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))
  
  ph <- ggplot(df_merged, aes(K,pkn,color =Process_type))+geom_bar(aes(linetype=Type),size = 0.7, stat="identity",alpha =0.0, position = "identity", fill= "white")+
    geom_vline(xintercept=3,  linetype="dashed")+ylab("Density")+xlab(TeX('$K_n$'))+
  #  ggtitle(TeX(sprintf('Posterior distribution for the number of clusters for $\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig_df_mut$N[1],fig_df_mut$N[(max(fig_df_mut$K)+1)],fig_df_mut$N[(2*max(fig_df_mut$K)+1)],fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))+
    theme_minimal()+ scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$n$ = %3.f',fig_df_mut$N[1]),sprintf('$n$ = %3.f',fig_df_mut$N[(max(fig_df_mut$K)+1)]),sprintf('$n$ = %3.f',fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),sprintf('$n$ = %3.f',fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))))+
    scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    scale_linetype_manual(name = "Distribution",values = c(1, 2),guide = guide_legend(override.aes = list(linetype = c(1, 2),color = "black") ) )+
    facet_wrap(~Process_type,labeller = pkn_names)
  ph
  
  pdf(paste0(fig_path,"Figure4_alpha_r_",a_g,"_",b_g,"_2.pdf" ))
  plot(ph)
  dev.off()
  
  
  E_k_<- df_merged %>% group_by(Process_type,Type) %>% summarize(sum(pkn *c(1:K_)))
  write.table(E_k_, file = paste0(fig_path,"Figure4_alpha_r_2.csv"), sep = ",", col.names = NA,qmethod = "double")
  
  
  ## Posterior of the weights
  
  pkn.labs<- c("W_1","W_2","W_3","W_4")
  names(pkn.labs) <- c(paste0("n = ", fig_df_mut$N[1]),paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]),paste0("n = ", fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),paste0("n = ", fig_df_mut$N[(3*max(fig_df_mut$K)+1)]))
  pkn_names <- as_labeller(
    c(`W_1` = paste0("n = ", fig_df_mut$N[1]), `W_2` = paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]),`W_3` = paste0("n = ", fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),`W_4` = paste0("n = ", fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))
  
  weights_fig <-fig_df$weights 
  weights_fig$n <-as.factor(weights_fig$W_val) 
  
  weights_fig_thin<- weights_fig%>%  group_by(K,n, W_,W_val) %>%  filter(row_number() %% 5 == 1)
  pw <- ggplot(weights_fig_thin, aes(x=K, y=weights, group =K, fill = n )) + ylab("Weights")+xlab(TeX('$K_n$'))+scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    geom_boxplot(alpha=0.5) +scale_fill_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$n$ = %3.f',fig_df_mut$N[1]),sprintf('$n$ = %3.f',fig_df_mut$N[(max(fig_df_mut$K)+1)]),sprintf('$n$ = %3.f',fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),sprintf('$n$ = %3.f',fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))))+
  #  ggtitle(TeX(sprintf('Posterior distribution of the component weights,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig_df_mut$N[1],fig_df_mut$N[(max(fig_df_mut$K)+1)],fig_df_mut$N[(2*max(fig_df_mut$K)+1)],fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))+
    theme_minimal()+  ylim(0,0.65)+facet_wrap(~W_,labeller = pkn_names)
  pw
  
  pdf(paste0(fig_path,"Figure4_alpha_r_",a_g,"_",b_g,"_4.pdf" ))
  plot(pw)
  dev.off()
#  ns = 100
#  alpha_s<- rgamma(ns, a_g,b_g*K_)
#  alpha_s_mod<- replace(alpha_s, alpha_s< 10^(-10), 10^(-10)) #to avoid small values for alpha, which could lead to inf values in funcDP/PYmo
  df_alpha_prior = tibble(it = rep(1:length(fig_df$alpha[[1]]),4), alpha_post = rep(rgamma(length(fig_df$alpha[[1]]), a_g,b_g*K_),4), label= c(rep(1,length(fig_df$alpha[[1]])),rep(2,length(fig_df$alpha[[2]])),rep(3,length(fig_df$alpha[[3]])),rep(4,length(fig_df$alpha[[4]])))) 
  df_alpha_prior$type = rep("Prior", dim(df_alpha_prior)[1])
  R_alpha_list<- list()
  R_alpha_list <- lapply(fig_df$alpha[1:4], gelman_diag)
  df_alpha = tibble(it = rep(1:length(R_alpha_list),4), alpha_post = flatten_dbl(fig_df$alpha[1:4]), label= c(rep(1,length(fig_df$alpha[[1]])),rep(2,length(fig_df$alpha[[2]])),rep(3,length(fig_df$alpha[[3]])),rep(4,length(fig_df$alpha[[4]])))) 
  df_alpha$type = rep("Posterior", dim(df_alpha)[1])
  df_alpha = rbind(df_alpha_prior, df_alpha)
  
  
  p_alpha <- ggplot(df_alpha, aes(x=alpha_post, color = type)) + scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D")+ 
    geom_density()+facet_wrap(~label)+ xlim(0,2)+ theme_minimal()
  p_alpha
  # Add mean line
  pdf(paste0(fig_path,"Figure4_alpha_r_",a_g,"_",b_g,"_5.pdf" ))
  plot(p_alpha)
  dev.off()
}
