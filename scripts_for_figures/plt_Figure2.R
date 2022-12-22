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

#---------- B) Specification of the simulation and prior parameters -----------------------------------------------

plt_fig2<-function(input_file, c_v =c(0.1, 0.5, 1, 2) , alpha_list, n_list, fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure2/" ){

  fig_df <- loadRData(input_file)
  fig_df_mut <- fig_df$line%>%group_by(Process_type,N)%>%mutate(pkn =density/sum(density))
  K_= max(fig_df_mut$K)
  julia <- julia_setup()
  julia_library("GibbsTypePriors")
  julia_assign("al", alpha_list[1])
  E_k =  round(sum(round(julia_eval("Pkn_Dirichlet_mult.(1:10,20, 10, al)"),3)*c(1:10)),1)
  
  
   p_l = ggplot(fig_df_mut, aes(x=K, colour = fig_df_mut$Process_type)) +
    geom_line(aes(x=K, y = pkn))  +  ylab('')+
    ylab('') + ggtitle(TeX(sprintf('Posterior dist for the number of clusters for $\\N =(%3.f, %2.f, %2.f)$,$\\E_k =%.2f$, \\hat{R} =(%.3f,%.3f,%.3f,%.3f)',fig_df_mut$N[1],fig_df_mut$N[(max(fig_df_mut$K)+1)],fig_df_mut$N[(2*max(fig_df_mut$K)+1)],5,fig_df_mut$Rh[1],fig_df_mut$Rh[(max(fig_df_mut$K)+1)],fig_df_mut$Rh[(2*max(fig_df_mut$K)+1)],fig_df_mut$Rh[(2*max(fig_df_mut$K)+1)])))+
    theme_minimal() +scale_x_continuous(limits = c(1, max(fig_df_mut$K)), expand = c(0, 0),breaks= c(1,seq(0,max(fig_df_mut$K),length=5)))+
    theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10))+geom_vline(xintercept=3,  linetype="dashed")+
    scale_color_discrete(name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$N$=%3.f',fig_df_mut$N[1]),sprintf('$N$=%3.f',fig_df_mut$N[(max(fig_df_mut$K)+1)]),sprintf('$N$=%3.f',fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),sprintf('$N$=%3.f',fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))))
  p_l
  
  pdf(paste0(fig_path,"Figure2_E_k_",E_k,"_1.pdf" ))
  plot(p_l)
  dev.off()
  
  
  julia <- julia_setup()
  julia_library("GibbsTypePriors")
  julia_assign("al", alpha_list)
  julia_assign("Nl", n_list)
  julia_assign("K_bound", K_)
 
  df_prior = tibble(K= 1:10, 
                    Pkn_1 = round(julia_eval("Pkn_Dirichlet_mult.(1:K_bound,Nl[1], K_bound, al[1])"),3),
                    Pkn_2= round(julia_eval("Pkn_Dirichlet_mult.(1:K_bound,Nl[2], K_bound, al[2])"),3), 
                    Pkn_3 = round(julia_eval("Pkn_Dirichlet_mult.(1:K_bound,Nl[3], K_bound, al[3])"),3),
                    Pkn_4 = round(julia_eval("Pkn_Dirichlet_mult.(1:K_bound,Nl[4], K_bound,al[4])"),3))%>% gather(Process_type, pkn,Pkn_1:Pkn_4)
  
  df_prior$Type =rep("Prior", dim(df_prior)[1])
  fig_df_mut$Type= rep("Posterior", dim(fig_df_mut)[1])
  
  df_merged = rbind(df_prior,fig_df_mut[, c("K","Process_type", "pkn", "Type")] )
  
  pkn.labs<- c("Pkn_1","Pkn_2","Pkn_3","Pkn_4")
  names(pkn.labs) <- c(paste0("n = ", fig_df_mut$N[1]),paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]),paste0("n = ", fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),paste0("n = ", fig_df_mut$N[(3*max(fig_df_mut$K)+1)]))
  pkn_names <- as_labeller(
    c(`Pkn_1` = paste0("n = ", fig_df_mut$N[1]), `Pkn_2` = paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]),`Pkn_3` = paste0("n = ", fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),`Pkn_4` = paste0("n = ", fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))
  
  ph <- ggplot(df_merged, aes(K,pkn,color =Process_type))+geom_bar(aes(linetype=Type),size = 0.7, stat="identity",alpha =0.0, position = "identity", fill= "white")+
    geom_vline(xintercept=3,  linetype="dashed")+ylab("Density")+xlab(TeX('$K_n$'))+
 #   ggtitle(TeX(sprintf('Posterior distribution for the number of clusters for $\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig_df_mut$N[1],fig_df_mut$N[(max(fig_df_mut$K)+1)],fig_df_mut$N[(2*max(fig_df_mut$K)+1)],fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))+
    theme_minimal()+ scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$n$ = %3.f',fig_df_mut$N[1]),sprintf('$n$ = %3.f',fig_df_mut$N[(max(fig_df_mut$K)+1)]),sprintf('$n$ = %3.f',fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),sprintf('$n$ = %3.f',fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))))+
    scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    scale_linetype_manual(name = "Distribution",values = c(1, 2),guide = guide_legend(override.aes = list(linetype = c(1, 2),color = "black") ) )+
    facet_wrap(~Process_type,labeller = pkn_names)
  ph
  
  pdf(paste0(fig_path,"Figure2_E_k_",E_k,"_2.pdf" ))
  plot(ph)
  dev.off()
  
  
  E_k_<- df_merged %>% group_by(Process_type,Type) %>% summarize(sum(pkn *c(1:K_)))
  write.table(E_k_, file = paste0(fig_path,"Figure2_E_k_",E_k,"_2.csv"), sep = ",", col.names = NA,qmethod = "double")
  

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
  
  pdf(paste0(fig_path,"Figure2_E_k_",E_k,"_3.pdf" ))
  plot(pw)
  dev.off()
  
  ### MTM 
  
  c_vec =c_v
  
  df_k_20<- df_post(n_list[1], c_vec=c_vec,ind=1, post = fig_df  )
  df_k_200<- df_post(n_list[2], c_vec = c_vec,ind=2, post = fig_df  )
  df_k_2000<- df_post(n_list[3], c_vec=c_vec,ind=3, post = fig_df  )
  df_k_20000<- df_post(n_list[4],c_vec=c_vec,ind=4, post = fig_df  )
  df_fin<- rbind(df_k_20,df_k_200,df_k_2000,df_k_20000)
  
  
  df_fin<- rbind(df_k_20,df_k_200,df_k_2000,df_k_20000)
  pkn.labs<- n_list
  names(pkn.labs) <- c(paste0("n = ", fig_df_mut$N[1]),paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]),paste0("n = ", fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),paste0("n = ", fig_df_mut$N[(3*max(fig_df_mut$K)+1)]))
  pkn_names <- as_labeller(
    c(`20` = paste0("n = ", fig_df_mut$N[1]), `200` = paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]),`2000` = paste0("n = ", fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),`20000` = paste0("n = ", fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))
  
  
  df_fin_bar = df_fin %>%
    group_by(Process_type,N, pkn) %>%
    summarize(count=n())%>% mutate(pkn_dens =count/sum(count))
  K_ = max(fig_df_mut$K)
  
  
  pm <- ggplot(df_fin_bar, aes(pkn,pkn_dens,color =Process_type, linetype = Process_type))+ geom_bar(aes(linetype=Process_type),size = 0.7, stat="identity",alpha =0.0, position = "identity", fill= "white")+
    geom_vline(xintercept=3,  linetype="dashed")+ylab("Density")+xlab(TeX('$K_n$'))+
  #  ggtitle(TeX(sprintf('PD for the num of clusters for MTM,$\\c_vec =(%2.1f,%2.1f,%2.1f, %2.1f) $ ',c_vec[1],c_vec[2],c_vec[3],c_vec[4])))+
    scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D",name = TeX(sprintf('$c$')) ,labels=unname(TeX(c(sprintf('$c$=%3.1f', c_vec[1]),sprintf('$c$=%3.1f',c_vec[2]),sprintf('$c$=%3.f',c_vec[3]),sprintf('$c$=%3.f',c_vec[4])))))+
    theme_minimal()+ ylab("Density")+xlab(TeX('$\tilde{K}$')) +
    scale_linetype_manual(values=c("solid", "dashed","longdash","dotted"),name = TeX(sprintf('$c$')),labels=unname(TeX(c(sprintf('$c$=%3.1f', c_vec[1]),sprintf('$c$=%3.1f',c_vec[2]),sprintf('$c$=%3.f',c_vec[3]),sprintf('$c$=%3.f',c_vec[4]))))) +
    scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    facet_wrap(~N)
  pm
  
  pdf(paste0(fig_path,"Figure2_E_k_",E_k,"_4.pdf" ))
  plot(pm)
  dev.off()

}


#data =  loadRData("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_40000.RData")
#K_ = 10
#pk <- MCMC_function(data, e0=0.56/K_, K=K_, M=M_it, burnin=nburn)
#M_it = 15000
#nburn = 5000
