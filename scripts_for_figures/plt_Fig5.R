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
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")

#---------- B) Specification of the simulation and prior parameters -----------------------------------------------

#input_file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig5.RData"

plt_fig1<-function(input_file, c_vec =c(0.1, 0.5, 1, 2) , fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure5/" ){
  
  # fig_path = "../figures/Figure1/"
  fig_df <- loadRData(input_file)
  
  fig_df_mut <- fig_df$line%>%group_by(Process_type,N)%>%mutate(pkn =density/sum(density))
  
  K_ = max(fig_df_mut$K)
  n_vec = c( fig_df_mut$N[1],fig_df_mut$N[(max(fig_df_mut$K)+1)],fig_df_mut$N[(2*max(fig_df_mut$K)+1)],fig_df_mut$N[(3*max(fig_df_mut$K)+1)])
  alpha_py = fig_df_mut$alpha_py[1]
  sigma_py = fig_df_mut$sigma_py[1]
  
  p_l = ggplot(fig_df_mut, aes(x=K, colour = fig_df_mut$Process_type)) +
    geom_line(aes(x=K, y = pkn))  +  ylab('')+
    ylab('') + ggtitle(TeX(sprintf('Posterior dist for the number of clusters for $\\N =(%3.f, %3.f, %3.f,%3.f)$,$\\alpha =%.2f$,$\\sigma =%.2f$, \\hat{R} =(%.3f,%.3f,%.3f,%.3f)',n_vec[1],n_vec[2],n_vec[3],n_vec[4],alpha_py,sigma_py,fig_df_mut$Rh[1],fig_df_mut$Rh[(max(fig_df_mut$K)+1)],fig_df_mut$Rh[(2*max(fig_df_mut$K)+1)],fig_df_mut$Rh[(3*max(fig_df_mut$K)+1)])))+
    theme_minimal() +scale_x_continuous(limits = c(1, max(fig_df_mut$K)), expand = c(0, 0),breaks= c(1,seq(0,max(fig_df_mut$K),length=5)))+
    theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10))+geom_vline(xintercept=3,  linetype="dashed")+
    scale_color_discrete(name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$N$=%3.f',n_vec[1]),sprintf('$N$=%3.f',n_vec[2]),sprintf('$N$=%3.f',n_vec[3]),sprintf('$N$=%3.f',n_vec[4])))))
  p_l
  
  pdf(paste0(fig_path,"Figure5_alpha_",alpha_py,"_",sigma_py,"_1.pdf" ))
  plot(p_l)
  dev.off()
  
  julia <- julia_setup()
  julia_library("GibbsTypePriors")
  
  
  julia_assign("K_bound", K_)
  julia_assign("N1", n_vec[1])
  julia_assign("N2", n_vec[2])
#  julia_assign("N3", n_vec[3])
 # julia_assign("N4", n_vec[4])
  # alpha = e0*K
  julia_assign("alpha_py", fig_df$line$alpha_py[1])
  julia_assign("sigma_py", fig_df$line$sigma_py[1])
    #a = julia_eval("Pkn_Dirichlet_mult.(1:K_bound,N1, N1, alpha)")
  
  df_prior = tibble(K= 1:K_, 
                    Pkn_1 = round(julia_eval("Pkn_PYM.(1:K_bound,N1, K_bound, alpha_py,sigma_py)"),3),
                    Pkn_2= round(julia_eval("Pkn_PYM.(1:K_bound,N2, K_bound, alpha_py,sigma_py)"),3))%>% gather(Process_type, pkn,Pkn_1:Pkn_2)
            #        Pkn_3 = round(julia_eval("Pkn_PYM.(1:K_bound,N3, K_bound, alpha_py,sigma_py)"),3),
                   # Pkn_4 = round(julia_eval("Pkn_PYM.(1:K_bound,N4, K_bound, alpha_py,sigma_py)"),3))%>% gather(Process_type, pkn,Pkn_1:Pkn_4)
  
  df_prior$Type =rep("Prior", dim(df_prior)[1])
  fig_df_mut$Type= rep("Posterior", dim(fig_df_mut)[1])
  
  df_merged = rbind(df_prior,fig_df_mut[, c("K","Process_type", "pkn", "Type")] )
  
  #pkn.labs<- c("Pkn_1","Pkn_2","Pkn_3","Pkn_4")
  pkn.labs<- c("Pkn_1","Pkn_2")
 # names(pkn.labs) <- c(paste0("n = ", fig_df_mut$N[1]),paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]),paste0("n = ", fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),paste0("n = ", fig_df_mut$N[(3*max(fig_df_mut$K)+1)]))
#  pkn_names <- as_labeller(
 #   c(`Pkn_1` = paste0("n = ", fig_df_mut$N[1]), `Pkn_2` = paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]),`Pkn_3` = paste0("n = ", fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),`Pkn_4` = paste0("n = ", fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))
   names(pkn.labs) <- c(paste0("n = ", fig_df_mut$N[1]),paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]))
   pkn_names <- as_labeller(c(`Pkn_1` = paste0("n = ", fig_df_mut$N[1]), `Pkn_2` = paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1])))
  
  p_h<- ggplot(df_merged, aes(K,pkn,color =Process_type))+geom_bar(aes(linetype=Type),size = 0.7, stat="identity",alpha =0.0, position = "identity")+
    geom_vline(xintercept=3,  linetype="dashed")+ylab("Density")+xlab(TeX('$K_n$'))+
    ggtitle(TeX(sprintf('Posterior distribution for the number of clusters for $\\alpha =%.3f$, $\\sigma =%.3f$,$\\n =(%2.f,%2.f,%2.f, %2.f) $ ',alpha_py,sigma_py,fig_df_mut$N[1],fig_df_mut$N[(max(fig_df_mut$K)+1)],fig_df_mut$N[(2*max(fig_df_mut$K)+1)],fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))+
    theme_minimal()+ scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$n$=%3.f',fig_df_mut$N[1]),sprintf('$n$=%3.f',fig_df_mut$N[(max(fig_df_mut$K)+1)]),sprintf('$n$=%3.f',fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),sprintf('$n$=%3.f',fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))))+
    scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    scale_linetype_manual(name = "Distribution",values = c(1, 2),guide = guide_legend(override.aes = list(linetype = c(1, 2),color = "black") ) )+
    facet_wrap(~Process_type,labeller = pkn_names)
  p_h
  pdf(paste0(fig_path,"Figure5_alpha_",alpha_py,"_",sigma_py,"_2.pdf" ))
  plot(p_h)
  dev.off()
  

  
  E_k<- df_merged %>% group_by(Process_type,Type) %>% summarize(sum(pkn *c(1:K_)))
  write.table(E_k, file =paste0(fig_path,"Prior_Posterior_exp_alpha_",alpha_py,"_",sigma_py,"_2.csv"), sep = ",", col.names = NA,
              qmethod = "double")
  
  
  
  
  ## Posterior of the weights
  
  #pkn.labs<- c("W_1","W_2","W_3","W_4")
  #names(pkn.labs) <- c(paste0("n = ", n_vec[1]),paste0("n = ",n_vec[2]),paste0("n = ", n_vec[3]),paste0("n = ", n_vec[4]))
  #pkn_names <- as_labeller(
  #  c(`W_1` = paste0("n = ", n_vec[1]), `W_2` = paste0("n = ", n_vec[2]),`W_3` = paste0("n = ", n_vec[3]),`W_4` = paste0("n = ", n_vec[4])))
  
  
  pkn.labs<- c("W_1","W_2")
  names(pkn.labs) <- c(paste0("n = ", n_vec[1]),paste0("n = ",n_vec[2]))
  pkn_names <- as_labeller(
    c(`W_1` = paste0("n = ", n_vec[1]), `W_2` = paste0("n = ", n_vec[2])))
  
  
  weights_fig <-fig_df$weights 
  weights_fig$n <-as.factor(weights_fig$W_val) 
  
  weights_fig_thin<- weights_fig%>%  group_by(K,n, W_,W_val) %>%  filter(row_number() %% 5 == 1)
  
  
  pw <- ggplot(weights_fig_thin, aes(x=K, y=weights, group =K, fill = n )) + ylab("Weights")+xlab(TeX('$K_n$'))+scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    geom_boxplot(alpha=0.5) +scale_fill_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$n$ = %3.f',n_vec[1]),sprintf('$n$ = %3.f',n_vec[2]),sprintf('$n$ = %3.f',n_vec[3]),sprintf('$n$ = %3.f',n_vec[4])))))+
    ggtitle(TeX(sprintf('Posterior distribution of the component weights $\\alpha =%.3f$,$\\sigma =%.3f$,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',alpha_py,sigma_py,n_vec[1],n_vec[2],n_vec[3],n_vec[4])))+
    theme_minimal()+  ylim(0,0.65)+facet_wrap(~W_,labeller = pkn_names)
  pw
  
  pdf(paste0(fig_path,"Figure5_alpha_",alpha_py,"_",sigma_py,"_3.pdf" ))
  plot(pw)
  dev.off()
  
  ### MTM 
  
  df_k_20<- df_post(n_vec[1], c_vec=c_vec,ind=1, post = fig_df  )
  df_k_200<- df_post(n_vec[2], c_vec = c_vec,ind=2, post = fig_df  )
#  df_k_2000<- df_post(n_vec[3], c_vec=c_vec,ind=3, post = fig_df  )
#  df_k_20000<- df_post(n_vec[4],c_vec=c_vec,ind=4, post = fig_df  )
#  df_fin<- rbind(df_k_20,df_k_200,df_k_2000,df_k_20000)
  df_fin<- rbind(df_k_20,df_k_200)
  
 # pkn.labs<- n_vec
#  names(pkn.labs) <- c(paste0("n = ", n_vec[1]),paste0("n = ",  n_vec[2]),paste0("n = ",  n_vec[3]),paste0("n = ",  n_vec[4]))
 # pkn_names <- as_labeller(
#    c(`20` = paste0("n = ",  n_vec[1]), `200` = paste0("n = ",  n_vec[2]),`2000` = paste0("n = ",  n_vec[3]),`20000` = paste0("n = ",  n_vec[4])))
  
   pkn.labs<- n_vec[1:2]
   names(pkn.labs) <- c(paste0("n = ", n_vec[1]),paste0("n = ",  n_vec[2]))
   pkn_names <- as_labeller(c(`20` = paste0("n = ",  n_vec[1]), `200` = paste0("n = ",  n_vec[2])))
  
  
  df_fin_bar = df_fin %>%
    group_by(Process_type,N, pkn) %>%
    summarize(count=n())%>% mutate(pkn_dens =count/sum(count))
  K_ = max(fig_df_mut$K)
  alpha_ =alpha_py
  
  
  pm <- ggplot(df_fin_bar, aes(pkn,pkn_dens,color =Process_type, linetype = Process_type))+ geom_bar(aes(linetype=Process_type),size = 0.7, stat="identity",alpha =0.0, position = "identity", fill= "white")+
    geom_vline(xintercept=3,  linetype="dashed")+
    ggtitle(TeX(sprintf('PD for the num of clusters for MTM $\\alpha =%.3f$,$\\sigma =%.3f$,$\\c_vec =(%2.1f,%2.1f,%2.1f, %2.1f) $ ',alpha_py,sigma_py,c_vec[1],c_vec[2],c_vec[3],c_vec[4])))+
    scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D",name = TeX(sprintf('$c$')) ,labels=unname(TeX(c(sprintf('$c$=%3.1f', c_vec[1]),sprintf('$c$=%3.1f',c_vec[2]),sprintf('$c$=%3.f',c_vec[3]),sprintf('$c$=%3.f',c_vec[4])))))+
    theme_minimal()+ ylab("Density")+xlab(expression(tilde(K))) +
    scale_linetype_manual(values=c("solid", "dashed","longdash","dotted"),name = TeX(sprintf('$c$')),labels=unname(TeX(c(sprintf('$c$=%3.1f', c_vec[1]),sprintf('$c$=%3.1f',c_vec[2]),sprintf('$c$=%3.f',c_vec[3]),sprintf('$c$=%3.f',c_vec[4]))))) +
    scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    facet_wrap(~N,labeller = pkn_names)
  pm
  
  pdf(paste0(fig_path,"Figure1_alpha_",alpha_py,"_",sigma_py,"_4.pdf" ))
  plot(pm)
  dev.off()
}
