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

#---------- B) Specification of the simulation and prior parameters -----------------------------------------------

input_file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig15.RData"

simulatuion_function_DPM<- function(nu1=10,nu2= 10*10,ns,Sn,N_tr){
  alpha_s<- rgamma(ns, nu1,nu2)
  alpha_s_mod<- replace(alpha_s, alpha_s< 10^(-10), 10^(-10)) #to avoid small values for alpha, which could lead to inf values in funcDP/PYmo
  #sum_list<- sapply(alpha_s_mod,funct,n=Sn,N=N_tr)
  p_list<- matrix(NA, nrow =length(alpha_s_mod), ncol=10 )
  julia_library("GibbsTypePriors")
  for (i in 1:length(alpha_s_mod)){
    julia_assign("al", alpha_s_mod[i])
    julia_assign("n", Sn)
    p_list[i,] =  round(julia_eval("Pkn_Dirichlet_mult.(1:10,n, 10, al)"),3)
  }
 return(colMeans(p_list))
}




plt_fig2<-function(input_file, c_v =c(0.1, 0.5, 1, 2), n_list= c(20,200,2000,20000), fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure4/" ){
  
  # fig_path = "../figures/Figure1/"
  fig_df <- loadRData(input_file)
  fig_df_mut <- fig_df$line%>%group_by(Process_type,N)%>%mutate(pkn =density/sum(density))
  K_= max(fig_df_mut$K)
  julia <- julia_setup()
  julia_library("GibbsTypePriors")
  
  
  
  julia <- julia_setup()
  julia_library("GibbsTypePriors")
   
  df_prior = tibble(K= 1:10, 
                    Pkn_1 = simulatuion_function_DPM(ns = 100, Sn=20, N_tr = 10 ),
                    Pkn_2= simulatuion_function_DPM(ns = 100, Sn=200, N_tr = 10 ), 
                    Pkn_3 = simulatuion_function_DPM(ns = 100, Sn=2000, N_tr = 10 ),
                    Pkn_4 = simulatuion_function_DPM(ns = 100, Sn=20000, N_tr = 10 ))%>% gather(Process_type, pkn,Pkn_1:Pkn_4)
  
  df_prior$Type =rep("Prior", dim(df_prior)[1])
  fig_df_mut$Type= rep("Posterior", dim(fig_df_mut)[1])
  
  df_merged = rbind(df_prior,fig_df_mut[, c("K","Process_type", "pkn", "Type")] )
  
  pkn.labs<- c("Pkn_1","Pkn_2","Pkn_3","Pkn_4")
  names(pkn.labs) <- c(paste0("n = ", fig_df_mut$N[1]),paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]),paste0("n = ", fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),paste0("n = ", fig_df_mut$N[(3*max(fig_df_mut$K)+1)]))
  pkn_names <- as_labeller(
    c(`Pkn_1` = paste0("n = ", fig_df_mut$N[1]), `Pkn_2` = paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]),`Pkn_3` = paste0("n = ", fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),`Pkn_4` = paste0("n = ", fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))
  
  ph <- ggplot(df_merged, aes(K,pkn,color =Process_type))+geom_bar(aes(linetype=Type),size = 0.7, stat="identity",alpha =0.0, position = "identity", fill= "white")+
    geom_vline(xintercept=3,  linetype="dashed")+ylab("Density")+xlab(TeX('$K$'))+
    ggtitle(TeX(sprintf('Posterior distribution for the number of clusters for $\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig_df_mut$N[1],fig_df_mut$N[(max(fig_df_mut$K)+1)],fig_df_mut$N[(2*max(fig_df_mut$K)+1)],fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))+
    theme_minimal()+ scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$n$ = %3.f',fig_df_mut$N[1]),sprintf('$n$ = %3.f',fig_df_mut$N[(max(fig_df_mut$K)+1)]),sprintf('$n$ = %3.f',fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),sprintf('$n$ = %3.f',fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))))+
    scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    scale_linetype_manual(name = "Distribution",values = c(1, 2),guide = guide_legend(override.aes = list(linetype = c(1, 2),color = "black") ) )+
    facet_wrap(~Process_type,labeller = pkn_names)
  ph
  
  pdf(paste0(fig_path,"Figure4_alpha_r_2.pdf" ))
  plot(ph)
  dev.off()
  
  
  ph2 <- ggplot(df_merged, aes(K,pkn,color =Process_type))+geom_bar(aes(linetype=Type),size = 0.7, stat="identity",alpha =0.0, position = "identity", fill= "white")+
    geom_vline(xintercept=3,  linetype="dashed")+ylab("Density")+xlab(TeX('$K$'))+
    #ggtitle(TeX(sprintf('Posterior distribution for the number of clusters for $\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig_df_mut$N[1],fig_df_mut$N[(max(fig_df_mut$K)+1)],fig_df_mut$N[(2*max(fig_df_mut$K)+1)],fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))+
    theme_minimal()+ scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$n$ = %3.f',fig_df_mut$N[1]),sprintf('$n$ = %3.f',fig_df_mut$N[(max(fig_df_mut$K)+1)]),sprintf('$n$ = %3.f',fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),sprintf('$n$ = %3.f',fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))))+
    scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    scale_linetype_manual(name = "Distribution",values = c(1, 2),guide = guide_legend(override.aes = list(linetype = c(1, 2),color = "black") ) )+
    facet_wrap(~Process_type,labeller = pkn_names)
  ph2
  
  pdf(paste0(fig_path,"Figure4_alpha_r_2_.pdf" ))
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
    ggtitle(TeX(sprintf('Posterior distribution of the component weights,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig_df_mut$N[1],fig_df_mut$N[(max(fig_df_mut$K)+1)],fig_df_mut$N[(2*max(fig_df_mut$K)+1)],fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))+
    theme_minimal()+  ylim(0,0.65)+facet_wrap(~W_,labeller = pkn_names)
  pw
  
  pdf(paste0(fig_path,"Figure2_E_k_",E_k,"_3.pdf" ))
  plot(pw)
  dev.off()
  
  
  ## Weights sorted by angle
  
  W_p_sorted <- list()
  for (j in 1:4){
    Eta_sp<- matrix(NA, nrow =dim(fig_df$eta[[1]])[1],ncol =  dim(fig_df$eta[[1]])[2] )
    Mu_sp <- array(rep(1, dim(fig_df$mu[[1]])[1]* dim(fig_df$mu[[1]])[2]* dim(fig_df$mu[[1]])[3]), dim=c( dim(fig_df$mu[[1]])[1],  dim(fig_df$mu[[1]])[2],  dim(fig_df$mu[[1]])[3]))
    Sig_sp <- array(rep(1, dim(fig_df$sig[[1]])[1]* dim(fig_df$sig[[1]])[2]* dim(fig_df$sig[[1]])[3]*dim(fig_df$sig[[1]])[4]), dim=c(dim(fig_df$sig[[1]])[1],  dim(fig_df$sig[[1]])[2], dim(fig_df$sig[[1]])[3], dim(fig_df$sig[[1]])[4]))
    for (i in 1: dim(fig_df$eta[[j]])[1]){
      pts =  fig_df$mu[[j]][i,,]
      pts_angle = apply(pts, 2, angle_dist)
      pts_sorted = sort(pts_angle[1,],index.return = TRUE)
      Eta_sp[i,] = fig_df$eta[[j]][i,pts_sorted$ix]
      Sig_sp[i,,,] = fig_df$sig[[j]][1,,,pts_sorted$ix]
    }  
    W_p_sorted[[j]] = Eta_sp
  }
  n_l = n_list
  dfW = tibble(K= 1:K_)
  for (j in 1:4){
    name_ <- paste("W_sp", j, sep = "")
    dfW[,name_]<- rep(n_l[j],K_)
  }
  dfW_ = dfW%>% gather(W_, W_val,  paste("W_sp", 1, sep = ""):paste("W_sp", length(n_l), sep = ""))
  
  W_df_sorted <- do.call(cbind, W_p_sorted)
  dfW_post <- cbind(dfW_,t(W_df_sorted))
  
  dfW_post_ <- gather(dfW_post, key = "it",value="weights", 4: dim(dfW_post)[2])
  
  
  
  pkn.labs<- c("W_sp1","W_sp2","W_sp3","W_sp4")
  names(pkn.labs) <- c(paste0("n = ", fig_df_mut$N[1]),paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]),paste0("n = ", fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),paste0("n = ", fig_df_mut$N[(3*max(fig_df_mut$K)+1)]))
  pkn_names <- as_labeller(
    c(`W_sp1` = paste0("n = ", fig_df_mut$N[1]), `W_sp2` = paste0("n = ", fig_df_mut$N[max(fig_df_mut$K)+1]),`W_sp3` = paste0("n = ", fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),`W_sp4` = paste0("n = ", fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))
  
  weights_fig_sort <-dfW_post_
  weights_fig_sort$n <-as.factor(weights_fig_sort$W_val) 
  
  weights_fig_sort_thin<- weights_fig_sort%>%  group_by(K,n, W_,W_val) %>%  filter(row_number() %% 5 == 1)
  
  
  pw2 <- ggplot(weights_fig_sort_thin, aes(x=K, y=weights, group =K, fill = n )) + ylab("Weights")+xlab(TeX('$K_n$'))+scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    geom_boxplot(alpha=0.5) +scale_fill_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$n$=%3.f',fig_df_mut$N[1]),sprintf('$n$=%3.f',fig_df_mut$N[(max(fig_df_mut$K)+1)]),sprintf('$n$=%3.f',fig_df_mut$N[(2*max(fig_df_mut$K)+1)]),sprintf('$n$=%3.f',fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))))+
    ggtitle(TeX(sprintf('Posterior distribution of the component weights ,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig_df_mut$N[1],fig_df_mut$N[(max(fig_df_mut$K)+1)],fig_df_mut$N[(2*max(fig_df_mut$K)+1)],fig_df_mut$N[(3*max(fig_df_mut$K)+1)])))+
    theme_minimal()+  facet_wrap(~W_,labeller = pkn_names)
  pw2
  
  pdf(paste0(fig_path,"Figure2_E_k_",E_k,"_3_2.pdf" ))
  plot(pw2)
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
    ggtitle(TeX(sprintf('PD for the num of clusters for MTM,$\\c_vec =(%2.1f,%2.1f,%2.1f, %2.1f) $ ',c_vec[1],c_vec[2],c_vec[3],c_vec[4])))+
    scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D",name = TeX(sprintf('$c$')) ,labels=unname(TeX(c(sprintf('$c$=%3.1f', c_vec[1]),sprintf('$c$=%3.1f',c_vec[2]),sprintf('$c$=%3.f',c_vec[3]),sprintf('$c$=%3.f',c_vec[4])))))+
    theme_minimal()+ ylab("Density")+xlab(TeX('$\tilde{K}$')) +
    scale_linetype_manual(values=c("solid", "dashed","longdash","dotted"),name = TeX(sprintf('$c$')),labels=unname(TeX(c(sprintf('$c$=%3.1f', c_vec[1]),sprintf('$c$=%3.1f',c_vec[2]),sprintf('$c$=%3.f',c_vec[3]),sprintf('$c$=%3.f',c_vec[4]))))) +
    scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    facet_wrap(~N)
  pm
  
  pdf(paste0(fig_path,"Figure2_E_k_",E_k,"_4.pdf" ))
  plot(pm)
  dev.off()
  
  n_vec =n_list
  
  df_MTM_MAP= tibble(c = c_vec,
                     Pkn_n1 = as.numeric(colnames(table(df_k_20$Process_type, df_k_20$pkn))[apply(table(df_k_20$Process_type, df_k_20$pkn), 1, which.max)]),
                     Pkn_n2=  as.numeric(colnames(table(df_k_200$Process_type, df_k_200$pkn))[apply(table(df_k_200$Process_type, df_k_200$pkn), 1, which.max)]),
                     Pkn_n3=  as.numeric(colnames(table(df_k_2000$Process_type, df_k_2000$pkn))[apply(table(df_k_2000$Process_type, df_k_2000$pkn), 1, which.max)]),
                     Pkn_n4=  as.numeric(colnames(table(df_k_20000$Process_type, df_k_20000$pkn))[apply(table(df_k_20000$Process_type, df_k_20000$pkn), 1, which.max)]))%>% gather(Process_type, pkn,Pkn_n1:Pkn_n4)
  
  
  df_MTM_MAP$name<- rep("MAP", dim(df_MTM_MAP)[1])
  df_MTM_Mean= tibble(c = c_vec,
                      Pkn_n1 = aggregate( df_k_20$pkn, list(df_k_20$Process_type), FUN=mean)$x,
                      Pkn_n2= aggregate( df_k_200$pkn, list(df_k_20$Process_type), FUN=mean)$x,
                      Pkn_n3=  aggregate( df_k_2000$pkn, list(df_k_20$Process_type), FUN=mean)$x,
                      Pkn_n4=  aggregate( df_k_20000$pkn, list(df_k_20$Process_type), FUN=mean)$x)%>% gather(Process_type, pkn,Pkn_n1:Pkn_n4)
  
  
  df_MTM_Mean$name<- rep("Mean", dim(df_MTM_MAP)[1])
  
  df_MTM <- rbind(df_MTM_MAP, df_MTM_Mean)
  
  
  write.table(df_MTM, file =paste0(fig_path,"Posterior_MTM_E_k_",E_k,".csv"), sep = ",", col.names = NA,
              qmethod = "double")
}
