
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

df_post_2<- function(c_vec, post , ind,  alpha, n=20){
  df_post_k= tibble(с= c_vec)
  df_post_k$val <- rep(NA, length(c_vec))
  for (i in 1:length(c_vec)){
    df_post_k$val[i] =  mean(round(apply_MTM(post$eta[[ind]],post$mu[[ind]],post$sig[[ind]], c =c_vec[i], n=n),3))
  }
  #df_mut = df_post_k%>%gather(Process_type, pkn,2:21)
  df_post_k$n =rep(n, length(c_vec))
  df_post_k$alpha= rep(alpha, length(c_vec))
  return(df_post_k)
}

df_post_3<- function(c_vec, post , ind,  alpha, n=20){
  df_post_k= tibble(с= c_vec)
  df_post_k$val <- rep(NA, length(c_vec))
  for (i in 1:length(c_vec)){
    result_MTM = apply_MTM(post$eta[[ind]],post$mu[[ind]],post$sig[[ind]], c =c_vec[i], n=n)
    map = names(table(result_MTM)[which.max(table(result_MTM))])
    df_post_k$val[i] = as.numeric(map)
  }
  #df_mut = df_post_k%>%gather(Process_type, pkn,2:21)
  df_post_k$n =rep(n, length(c_vec))
  df_post_k$alpha= rep(alpha, length(c_vec))
  return(df_post_k)
}



plt_fig3<-function(input_file, c_vec,alpha=0.01, fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure3/" ){
   fig_df <- loadRData(input_file)
   fig_df_mut <- fig_df$line%>%group_by(Process_type,N)%>%mutate(pkn =density/sum(density))
   K_ = max(fig_df_mut$K)
   n_vec  = c( fig_df_mut$N[1],fig_df_mut$N[max(fig_df_mut$K)+1], fig_df_mut$N[(2*max(fig_df_mut$K)+1)],fig_df_mut$N[(3*max(fig_df_mut$K)+1)])
  ### MTM 
  df_1<-  df_post_2(c_vec, fig_df, ind = 1, alpha = alpha , n=n_vec[1])
  df_2<- df_post_2(c_vec, fig_df,  ind = 2, alpha = alpha , n=n_vec[2])
  df_3<-  df_post_2(c_vec, fig_df,  ind = 3, alpha = alpha , n=n_vec[3])
  df_4<- df_post_2(c_vec,fig_df,  ind = 4, alpha = alpha , n=n_vec[4])

  df_fin= rbind(df_1, df_2, df_3, df_4)
  df_fin$type = rep("Mean",dim(df_fin)[1])
  
  df_1_map<-  df_post_3(c_vec, fig_df, ind = 1, alpha = alpha , n=n_vec[1])
  df_2_map<- df_post_3(c_vec, fig_df,  ind = 2, alpha = alpha , n=n_vec[2])
  df_3_map<-  df_post_3(c_vec, fig_df,  ind = 3, alpha = alpha , n=n_vec[3])
  df_4_map<- df_post_3(c_vec,fig_df,  ind = 4, alpha = alpha , n=n_vec[4])

  df_fin_map =  rbind(df_1_map, df_2_map, df_3_map, df_4_map)
  df_fin_map$type = rep("MAP",dim(df_fin)[1])
  
  df_fin<- rbind( df_fin_map,df_fin)
pkn.labs<- n_vec
names(pkn.labs) <- c(paste0("n = ", n_vec[1]),paste0("n = ", n_vec[2]),paste0("n = ", n_vec[3]),paste0("n = ",  n_vec[4]))
pkn_names <- as_labeller(
  c(`20` = paste0("n = ",  n_vec[1]), `200` = paste0("n = ",  n_vec[2]),`2000` = paste0("n = ",  n_vec[3]),`20000` = paste0("n = ",  n_vec[4])))
df_fin$n = as.factor(df_fin$n)

pm <- ggplot(df_fin, aes(x = с, y = val, color = n))+ geom_line(aes(linetype=type))+
  theme_minimal()+ ylab(expression(tilde(K)))+xlab(TeX('$c$')) +geom_hline(yintercept = 3, linetype="dashed", size = 0.5, alpha =0.5) +
  scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D",name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$n$ = %3.f', n_vec[1]),sprintf('$n$ = %3.f',n_vec[2]),sprintf('$n$ = %3.f',n_vec[3]),sprintf('$n$ = %3.f',n_vec[4])))))+
  scale_linetype_manual(name = "Distribution",values = c(4, 1),guide = guide_legend(override.aes = list(linetype = c(4, 1),color = "black") ) )+
  facet_wrap(~n,labeller = pkn_names)
pm
pdf(paste0(fig_path,"Figure3_alpha_",alpha,"_6.pdf" ))
plot(pm)
dev.off()
return(df_fin)
}
