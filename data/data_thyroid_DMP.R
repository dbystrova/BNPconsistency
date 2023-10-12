rm(list=ls()) 
setwd("~/Documents/GitHub/BNPconsistency/data")

library(tidyverse)
library(latex2exp)
library(mclust)

data("thyroid")
save(thyroid, file = "~/Documents/GitHub/BNPconsistency/data/thyroid.RData")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Random_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Estimation_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Identification_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Gibbs_sampling_function.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/plt_Figure1.R")

data_raw = loadRData("~/Documents/GitHub/BNPconsistency/data/thyroid.RData")
data<- list()
data$y = data_raw[,2:6]

comparison_data_alpha<- function(data,K_, M_it, nburn, alpha_l){
  pk<- list()
  R_h <- c()
  W_non_sorted <- list()
  W <- list()
  Mu_mat <- list()
  S_mat<- list()
 
  N<- dim(data$y)[1]
  for (i in 1:length(alpha_l)){
    pk[[i]] <- MCMC_function(data, e0=alpha_l[i], K=K_, M=M_it, burnin=nburn)  #conversion from Dir to DPM 
    Eta_<- matrix(NA, nrow =dim(pk[[i]]$Eta)[1],ncol =  dim(pk[[i]]$Eta)[2] )
    for (j in 1:dim(pk[[i]]$Eta)[1]){ Eta_[j,] <- sort(pk[[i]]$Eta[j,],decreasing = TRUE)}
    W[[i]] <- Eta_
    W_non_sorted[[i]] <- pk[[i]]$Eta
    Mu_mat[[i]] <- pk[[i]]$Mu
    S_mat[[i]] <-pk[[i]]$Sigma
  }
  df_hist = tibble(K= 1:K_)
  df_info = tibble(models= 1:length(pk))
  df_weights = tibble(K= 1:K_)
  df_pk_all = tibble(it= 1:((M_it -nburn)*2))
  for (j in 1:length(alpha_l)){
    name_ <- paste("Pkn_", j, sep = "")
    df_info[j,"Rh"] =pk[[j]]$ll_rhat
    df_info[j,"N"] =N
    name5_ <- paste("Alpha_", j, sep = "")
    df_hist[,name_]<- pk[[j]]$p_k
    df_weights[,name5_]<- rep(alpha_l[j],length(pk[[j]]$p_k))
    df_pk_all[,name_]<- pk[[j]]$p_k_all
  }
  df = df_hist%>% gather(Process_type, density,  paste("Pkn_", 1, sep = ""):paste("Pkn_", length(alpha_l), sep = ""))
  
  df5 = df_weights%>% gather(Alpha, Alpha_val,  paste("Alpha_", 1, sep = ""):paste("Alpha_", length(alpha_l), sep = ""))
  W_df <- do.call(cbind, W)
  df5_post <- cbind(df5,t(W_df))
  df5_post_ <- gather(df5_post, key = "it",value="weights", 4: dim(df5_post)[2])
  
  df_pk_all_ = df_pk_all%>% gather(Process_type, density,  paste("Pkn_", 1, sep = ""):paste("Pkn_", length(alpha_l), sep = ""))
  
    write.table(df_info, file =paste0(fig_path,paste0("Convergence_res_data_it_",M_it,"burn_",nburn,".csv")), sep = ",", col.names = NA,
              qmethod = "double")
  return(list(line = df, hist =df_pk_all_, weights = df5_post_, eta = W_non_sorted, mu = Mu_mat, sig = S_mat))
}


data_file="~/Documents/GitHub/BNPconsistency/data/thyroid.RData"
data_raw =  loadRData(data_file)
data<- list()
data$y = data_raw[,2:6]
N = dim(data$y)[1]
alpha_l = c(0.5, 1, 2, 8)
fig_path=''
final = comparison_data_alpha(data= data,K_= 10 , M_it = 20000, nburn= 10000, alpha_l=alpha_l)
save(final, file = "~/Documents/GitHub/BNPconsistency/data/df_thyroid.RData")

#final = loadRData( "~/Documents/GitHub/BNPconsistency/data/df_thyroid2.RData")

## Normal plots

 # fig_path = "../figures/Figure1/"
  fig_df <- final
  
  fig_df_mut <- fig_df$line%>%group_by(Process_type)%>%mutate(pkn =density/sum(density))
  
  
  julia <- julia_setup()
  julia_library("GibbsTypePriors")
  
  K_= 10
  julia_assign("K_bound", K_)
  julia_assign("A1", alpha_l[1]*K_)
  julia_assign("A2", alpha_l[2]*K_)
   julia_assign("A3", alpha_l[3]*K_)
   julia_assign("A4", alpha_l[4]*K_)
  # alpha = e0*K
  julia_assign("N", N)
 # print( fig_df$line$alpha[1]*K_)
  #a = julia_eval("Pkn_Dirichlet_mult.(1:K_bound,N,K_bound, A1)")
  
  df_prior = tibble(K= 1:K_, 
                    Pkn_1 = round(julia_eval("Pkn_Dirichlet_mult.(1:10,215, 10, A1)"),3),
                    Pkn_2= round(julia_eval("Pkn_Dirichlet_mult.(1:10,215, 10, A2)"),3),
                     Pkn_3 = round(julia_eval("Pkn_Dirichlet_mult.(1:10,215, 10, A3)"),3),
                     Pkn_4 = round(julia_eval("Pkn_Dirichlet_mult.(1:10,215, 10, A4)"),3)
                    )%>% gather(Process_type, pkn,Pkn_1:Pkn_4)
  
  df_prior$Type =rep("Prior", dim(df_prior)[1])
  fig_df_mut$Type= rep("Posterior", dim(fig_df_mut)[1])
  
  df_merged = rbind(df_prior,fig_df_mut[, c("K","Process_type", "pkn", "Type")] )
  
  
  
  pkn.labs<- c("Pkn_1","Pkn_2","Pkn_3","Pkn_4")
  names(pkn.labs) <- c(paste0("alpha = ", alpha_l[1]),paste0("alpha= ", alpha_l[1]))
  pkn_names <- as_labeller(c(`Pkn_1` = paste0(TeX("alpha = "), alpha_l[1]), `Pkn_2` = paste0("alpha = ", alpha_l[2]),`Pkn_3` = paste0("alpha = ", alpha_l[3]),`Pkn_4` = paste0("alpha = ", alpha_l[4])))
  
  
  
  p_h<- ggplot(df_merged, aes(K,pkn,color =Process_type))+geom_bar(aes(linetype=Type),size = 0.7, stat="identity",alpha =0.0, position = "identity")+
    geom_vline(xintercept=3,  linetype="dashed")+ylab("Density")+xlab(TeX('$K_n$'))+
  #  ggtitle(TeX(sprintf('Posterior distribution for the number of clusters')))+
    theme_minimal()+ scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$alpha$')))+
    scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    scale_linetype_manual(name = "Distribution",values = c(1, 2),guide = guide_legend(override.aes = list(linetype = c(1, 2),color = "black") ) )+
    facet_wrap(~Process_type,labeller = pkn_names)
  p_h
  pdf(paste0(fig_path,"Figure_data_1.pdf" ))
  plot(p_h)
  dev.off()
  
   E_k<- df_merged %>% group_by(Process_type,Type) %>% summarize(sum(pkn *c(1:K_)))
  write.table(E_k, file =paste0(fig_path,"Prior_Posterior_exp_alpha_data.csv"), sep = ",", col.names = NA,
              qmethod = "double")
  
  
  
  
  ## Posterior of the weights
  
 
  weights_fig <-final$weights
  weights_fig$Alpha<-as.factor(weights_fig$Alpha) 
  
  weights_fig_thin<- weights_fig%>%  group_by(K, Alpha,Alpha_val) %>%  filter(row_number() %% 5 == 1)
  
  
  pw <- ggplot(weights_fig_thin, aes(x=K, y=weights, group =K )) + ylab("Weights")+xlab(TeX('$K_n$'))+scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    geom_boxplot(alpha=0.5) +scale_fill_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D")+
   # ggtitle(TeX(sprintf('Posterior distribution of the component weights $\\alpha =%.3f$,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig_df_mut$alpha[1],n_vec[1],n_vec[2],n_vec[3],n_vec[4])))+
    theme_minimal()+  ylim(0,0.65)+facet_wrap(~Alpha_val)
  pw
  
  pdf(paste0(fig_path,"Figure_data__2.pdf" ))
  plot(pw)
  dev.off()
 

####MTM 

 # function(p_eta, p_mu, p_sig,c,n, w_n = (log(n)/n)^(1/4))

df_post_MTM_alpha<- function(c_vec, post,  alpha_l, N){
  df_post_k= tibble(с= c_vec)
  n_it = dim( post$eta[[1]])[1]
  tokeep <- seq(2, n_it, 4)
  for  (i in 1:length(alpha_l)){
    post$eta[[i]] =   post$eta[[i]][tokeep,]
    post$mu[[i]] = post$mu[[i]][tokeep,,] 
    post$sig[[i]] = post$sig[[i]][tokeep,,,] 
  }
  for (i in 1:length(c_vec)){
    for (j in 1:length(alpha_l)){
      df_post_k[ i,(j+1)] = mean(round(apply_MTM(p_eta = post$eta[[j]],p_mu = post$mu[[j]],p_sig = post$sig[[j]], c =c_vec[i], n=N),3))
    }
  }
  names(df_post_k) = c('c', as.character(alpha_l))
  df_mut_mean = df_post_k%>%gather(key =alpha_level,value =val, 2:(length(alpha_l)+1))
  df_post_k_map= tibble(с= c_vec)
  df_mut_mean$type = rep("Mean", dim(df_mut_mean)[1])
  for (i in 1:length(c_vec)){
    for (j in 1:length(alpha_l)){
      result_MTM = apply_MTM(post$eta[[j]],post$mu[[j]],post$sig[[j]], c =c_vec[i], n=N)
      map = names(table(result_MTM)[which.max(table(result_MTM))])
      df_post_k_map[i,(j+1)]  = as.numeric(map)
     }
  }
  names(df_post_k_map) = c('c', as.character(alpha_l))
  df_mut_map = df_post_k_map%>%gather(key =alpha_level,value =val, 2:(length(alpha_l)+1))
  df_mut_map$type = rep("MAP", dim(df_mut_map)[1])
  df_mut = rbind(df_mut_mean,df_mut_map)
  
  
   df_mut$alpha_level = as.factor(df_mut$alpha_level)
  
  pm <- ggplot(df_mut, aes(x = c, y = val, color = alpha_level))+ geom_line(aes(linetype=type))+
    theme_minimal()+ ylab(expression(tilde(K)))+xlab(TeX('$c$')) +geom_hline(yintercept = 3, linetype="dashed", size = 0.5, alpha =0.5) +    facet_wrap(~alpha_level)
  plot(pm)
  pdf(paste0("Figure_data_MtM.pdf"))
  plot(pm)
  dev.off()
}

c_vec_long =seq(0, 1, length.out = 15)

df_post_MTM_alpha(c_vec = c_vec_long, post= final,  alpha_l, N = 215)
