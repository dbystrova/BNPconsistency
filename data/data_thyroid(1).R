rm(list=ls()) 
setwd("~/Documents/GitHub/BNPconsistency/data")

#### Try on BNPmix

library(BNPmix)
library(tidyverse)
library(latex2exp)
library(viridis)
library(JuliaCall)
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")

library(BSDA)
library(mclust)
data("thyroid")
data("Slc")


l2norm_MV_gaussian = function(theta.j, theta.i, diag = FALSE){
  l2_mu = sum((theta.j[[1]] -theta.i[[1]] )^2)
  #l2_Sigma = sqrt(sum((theta.j[[2]] -theta.i[[2]] )^2))
  #l2_Sigma =norm(theta.j[[2]] -theta.i[[2]],"2")^2
  if (diag){  l2_Sigma = sum((theta.j[[2]] -theta.i[[2]] )^2) }
  else{
  # l2_Sigma = sum(((theta.j[[2]] - theta.i[[2]])[upper.tri(theta.j[[2]] -theta.i[[2]], diag = TRUE)])^2)}
   l2_Sigma =norm(theta.j[[2]] - theta.i[[2]],"F")^2}
  return(sqrt(l2_mu + l2_Sigma ))
}

l2norm_univ = function(theta.j, theta.i){
  l2_mu = sum((theta.j[[1]] -theta.i[[1]])^2)
  l2_Sigma = sum((theta.j[[2]] -theta.i[[2]])^2)
  return(sqrt(l2_mu + l2_Sigma))
}

apply_MTM_univ <- function(p_eta, p_mu, p_sig,c,n, w_n = (log(n)/n)^(1/4)){
  post_k<-c()
  it <- dim(p_eta)[1]
  for (k in 1:it){
    G <- list()
    G$p = p_eta[k,]
    G$theta = list ()
    for (i in 1:length(p_eta[k,])){
      G$theta[[i]]<-  list(p_mu[k,i],p_sig[k,i])
    }
    G.post <- MTM_univ(G, w_n, c)
    #print(G.post$p)
    post_k[k]<- length(G.post$p)
  }
  
  return(post_k)
}

MTM_univ <- function(G, w_n, c){
  p.k = G$p
  theta.k= G$theta
  ## stage 1
  #print(p.k)
  if (length( which(p.k == 0))>0) {
    theta.k = theta.k[-c(which(p.k == 0))]
    p.k =p.k[-c(which(p.k == 0))]
  }
  theta_ind = sample(1:length(theta.k), length(theta.k), replace=FALSE, prob = p.k)
  tau = theta_ind
  #print(tau)
  p_new = vector(mode="numeric", length=length(p.k ))
  p_new= p.k
  t = length(tau)
  i = 1 
  while (i<=t){
    #print(tau[i])
    j = 1
    while (j<=t){
      #print(tau[j])
      if(tau[j]< tau[i]){
        if(l2norm_univ(theta.k[[tau[j]]], theta.k[[tau[i]]]) <= w_n){
          p_new[tau[i]] =p.k[tau[j]] + p.k[tau[i]]
          # p_new = p_new[-tau[j]]
          tau = tau[-j]
          t = length(tau)
          #print(c(i,j))
          if (j<i) {i = i-1}
          j= j-1
          #print(c(i,j))
        }
      }
      j = j + 1
    }
    i = i +1
  }

  # cat("tau : ", tau, "length = ", length(tau), "\n")
  #stage 2
  sorted_weights  = sort(p_new[tau], decreasing = TRUE, index.return=TRUE)
  # cat("sort_ind : ", sorted_weights$ix, "length = ", length(sorted_weights$ix), "\n")
  p = sorted_weights$x
  th = theta.k[tau]
  th_sorted = th[sorted_weights$ix]
  # p = p_sorted[tau]
  
  A = which(p > (c*w_n)^2)
  N =  which(p <= (c*w_n)^2)
  
  for(i in A){
    for(j in A){
      if (j < i){
        #  print(p[i]*(l2norm_MV_gaussian(th[[i]],th[[j]]))^2)
        if (p[i]*(l2norm_univ(th[[i]],th[[j]]))^2 <= (c*w_n)^2){
          N = c(N, i)
          A = A[-i]
          #    print(A)
        }
      }
    }
  }
  
  for(i in N){
    v=0
    ind_min = 1
    for (j in A){
      if (v<= l2norm_univ(th[[j]],th[[i]])){
        v = l2norm_univ(th[[j]],th[[i]])
        ind_min = j
      }
    }
    p[ind_min] = p[ind_min] + p[i] 
  }
  return (list(p=p[A], th=th[A]))
}

###### DMP
save(thyroid, file = "~/Documents/GitHub/BNPconsistency/data/thyroid.RData")
save(Slc, file = "~/Documents/GitHub/BNPconsistency/data/Slc.RData")

source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Random_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Estimation_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Identification_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Gibbs_sampling_function.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/plt_Figure1.R")


df_post_MTM_alpha<- function(c_vec, post,  alpha_l, N, name){
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
      df_post_k[ i,(j+1)] = mean(round(apply_MTM(post$eta[[j]],post$mu[[j]],post$sig[[j]], c =c_vec[i], n=N),3))
    }
  }
  names(df_post_k) = c('c', as.character(alpha_l))
  df_mut_mean = df_post_k%>%gather(key =alpha_level,value =val, 2:5)
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
  df_mut_map = df_post_k_map%>%gather(key =alpha_level,value =val, 2:5)
  df_mut_map$type = rep("MAP", dim(df_mut_map)[1])
  df_mut = rbind(df_mut_mean,df_mut_map)
  
  
  df_mut$alpha_level = as.factor(df_mut$alpha_level)
  
  pm <- ggplot(df_mut, aes(x = c, y = val, color = alpha_level))+ geom_line(aes(linetype=type))+
    scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D",name = TeX(sprintf('$\\bar{\\alpha}$')), labels=unname(TeX(c(sprintf('$\\bar{\\alpha}$ = %.2f', alpha_l[1]),
                                                                                                                                  sprintf('$\\bar{\\alpha}$ = %.1f', alpha_l[2]),
                                                                                                                                  sprintf('$\\bar{\\alpha}$ = %.1f', alpha_l[3]),
                                                                                                                                  sprintf('$\\bar{\\alpha}$ = %1.f', alpha_l[4])))))+
    theme_minimal()+theme(strip.background = element_blank(), strip.text.x = element_blank())+
    ylab(expression(tilde(K)))+xlab(TeX('$c$')) +geom_hline(yintercept = 3, linetype="dashed", size = 0.5, alpha =0.5) +    facet_wrap(~alpha_level)
  plot(pm)
  pdf(paste0("Figure_thyroid_DPM_MTM.pdf"))
  plot(pm)
  dev.off()
}

df_post_univ <- function(c_vec, post , ind,  alpha, N){
  df_post_k= tibble(с= c_vec)
  df_post_k$val <- rep(NA, length(c_vec))
  df_post_k_map= tibble(с= c_vec)
  df_post_k_map$val <- rep(NA, length(c_vec))
  
  n_it = dim( post$eta[[1]])[1]
  tokeep <- seq(2, n_it, 4)
  post$eta[[ind]] =   post$eta[[ind]][tokeep,]
  post$mu[[ind]] = post$mu[[ind]][tokeep,,] 
  post$sig[[ind]] = post$sig[[ind]][tokeep,,,] 
  
  for (i in 1:length(c_vec)){
    result_MTM = apply_MTM_univ(post$eta[[ind]],post$mu[[ind]],post$sig[[ind]], c = c_vec[i], n=N)
    map = names(table(result_MTM)[which.max(table(result_MTM))])
    df_post_k_map$val[i] = as.numeric(map)
    df_post_k$val[i] = mean(round(result_MTM,3))
  }
  #df_mut = df_post_k%>%gather(Process_type, pkn,2:21)
  df_post_k$alpha= rep(alpha, length(c_vec))
  df_post_k$type = rep("Mean", length(c_vec))
  df_post_k_map$alpha= rep(alpha, length(c_vec))
  df_post_k_map$type = rep("MAP", length(c_vec))
  return(list(df_post_k=df_post_k,df_post_k_map=df_post_k_map))
}


plt_univ_MTM<- function(c_list, post,  alpha_l, N){
  df_1 <- df_post_univ(c_list[[1]], post, ind = 1, alpha = alpha_l[1], N)
  df_2 <- df_post_univ(c_list[[2]], post, ind = 2, alpha = alpha_l[2], N)
  df_3 <- df_post_univ(c_list[[3]], post, ind = 3, alpha = alpha_l[3], N)
  df_4 <- df_post_univ(c_list[[4]], post, ind = 4, alpha = alpha_l[4], N)
  
  df_fin= rbind(df_1$df_post_k, df_2$df_post_k, df_3$df_post_k, df_4$df_post_k)
  df_fin_map= rbind(df_1$df_post_k_map, df_2$df_post_k_map, df_3$df_post_k_map, df_4$df_post_k_map)
  
  df_mut <- rbind(df_fin_map,df_fin)
  df_mut$alpha = as.factor(df_mut$alpha)
  
  # df_mut$alpha_level = as.factor(df_mut$alpha_level)
  return(df_mut)
}

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
  df_ = tibble(K= 1:K_)
  df4_ = tibble(K= 1:K_)
  df5_ = tibble(K= 1:K_)
  for (j in 1:length(alpha_l)){
    name_ <- paste("Pkn_", j, sep = "")
    name4_ <- paste("Alpha_", j, sep = "")
    name5_ <- paste("W_", j, sep = "")
    df_[,name_]<- pk[[j]]$p_k
    print(pk[[j]]$ll_rhat)
    df4_[,name4_]<- rep(alpha_l[j]/K_,length(pk[[j]]$p_k))
    df5_[,name5_]<- rep(N[j],length(pk[[j]]$p_k))
  }
  df = df_%>% gather(Process_type, density,  paste("Pkn_", 1, sep = ""):paste("Pkn_", length(alpha_l), sep = ""))
  df4 = df4_%>% gather(Alpha_, Alpha_val,  paste("Alpha_", 1, sep = ""):paste("Alpha_", length(alpha_l), sep = ""))
  df5 = df5_%>% gather(W_, W_val,  paste("W_", 1, sep = ""):paste("W_", length(alpha_l), sep = ""))
  
  W_df <- do.call(cbind, W)
  df5_post <- cbind(df5,t(W_df))
  
  df5_post_ <- gather(df5_post, key = "it",value="weights", 4: dim(df5_post)[2])
  
  df$Al =df4$Alpha_val
  
  df_l_ = tibble(K= 1:((M_it)*2))
  df4_l_ = tibble(K= 1:((M_it)*2))
  for (j in 1:length(alpha_l)){
    name_ <- paste("P_", j, sep = "")
    name4_ <- paste("Alpha_", j, sep = "")
    df_l_[,name_]<- pk[[j]]$p_k_all
    df4_l_[,name4_]<- rep(alpha_l[j],length(pk[[j]]$p_k_all))
  }
  df_l = df_l_%>% gather(Process_type, density,  paste("P_", 1, sep = ""):paste("P_", length(alpha_l), sep = ""))
  df_l4 = df4_l_%>% gather(Alpha_, Alpha_val,  paste("Alpha_", 1, sep = ""):paste("Alpha_", length(alpha_l), sep = ""))
  df_l$Al =df_l4$Alpha_val
  return(list(line = df, hist =df_l, weights = df5_post_, eta = W_non_sorted, mu = Mu_mat, sig = S_mat))
}


##################################
######## Multivariate DMP ########
##################################

data_raw = loadRData("~/Documents/GitHub/BNPconsistency/data/thyroid.RData")
data<- list()
data$y = data_raw[,2:6]
alpha = 1
K_ = 10
M_it = 50000
nburn = 30000
N =dim( data_raw[,2:6])[1]

alpha_l = c(0.01,0.5, 1, 10)
# final = comparison_data_alpha(data= data,K_= 10 , M_it = M_it, nburn= nburn, alpha_l=alpha_l)
# save(final, file = "~/Documents/GitHub/BNPconsistency/data/df_thyroid.RData")

final = loadRData( "~/Documents/GitHub/BNPconsistency/data/df_thyroid.RData")

## Normal plots

 # fig_path = "../figures/Figure1/"
fig_path = ""
fig_df <- final
  
  fig_df_mut <- fig_df$line%>%group_by(Process_type,Al)%>%mutate(pkn =density/sum(density))
  
  
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
                    Pkn_4 = round(julia_eval("Pkn_Dirichlet_mult.(1:10,215, 10, A4)"),3))%>% gather(Process_type, pkn,Pkn_1:Pkn_4)
  
  df_prior$Type =rep("Prior", dim(df_prior)[1])
  fig_df_mut$Type= rep("Posterior", dim(fig_df_mut)[1])
  
  df_merged = rbind(df_prior,fig_df_mut[, c("K","Process_type", "pkn", "Type")] )
  
  
  p_h<- ggplot(df_merged, aes(K,pkn,color =Process_type))+geom_bar(aes(linetype=Type),size = 0.7, stat="identity",alpha =0.0, position = "identity")+
    geom_vline(xintercept=3,  linetype="dashed")+ylab("Density")+xlab(TeX('$K_n$'))+
    theme_minimal()+theme(strip.background = element_blank(), strip.text.x = element_blank())+
    scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D",name = TeX(sprintf('$\\bar{\\alpha}$')), labels=unname(TeX(c(sprintf('$\\bar{\\alpha}$ = %.2f', alpha_l[1]),
                                                                                                                                                  sprintf('$\\bar{\\alpha}$ = %.1f', alpha_l[2]),
                                                                                                                                                  sprintf('$\\bar{\\alpha}$ = %.1f', alpha_l[3]),
                                                                                                                                                  sprintf('$\\bar{\\alpha}$ = %1.f', alpha_l[4])))))+
    scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
    scale_linetype_manual(name = "Distribution",values = c(1, 2),guide = guide_legend(override.aes = list(linetype = c(1, 2),color = "black") ) )+
    facet_wrap(~Process_type)
  p_h
  pdf(paste0(fig_path,"Figure1_DMP_thyroid.pdf" ))
  plot(p_h)
  dev.off()
  
   E_k<- df_merged %>% group_by(Process_type,Type) %>% summarize(sum(pkn *c(1:K_)))
  write.table(E_k, file =paste0(fig_path,"Prior_Posterior_exp_alpha_.csv"), sep = ",", col.names = NA,
              qmethod = "double")
  
  
  
  
  # ## Posterior of the weights
  # 
  # 
  # weights_fig <-v 
  # weights_fig$n <-as.factor(weights_fig$W_val) 
  # 
  # weights_fig_thin<- weights_fig%>%  group_by(K,n, W_,W_val) %>%  filter(row_number() %% 5 == 1)
  # 
  # 
  # pw <- ggplot(weights_fig_thin, aes(x=K, y=weights, group =K )) + ylab("Weights")+xlab(TeX('$K_n$'))+scale_x_continuous(breaks = c(1,3,7,10), limits = c(0,11))+
  #   geom_boxplot(alpha=0.5) +scale_fill_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D")+
  #  # ggtitle(TeX(sprintf('Posterior distribution of the component weights $\\bar{\\alpha} =%.3f$,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig_df_mut$alpha[1],n_vec[1],n_vec[2],n_vec[3],n_vec[4])))+
  #   theme_minimal()+  ylim(0,0.65)+facet_wrap(~W_)
  # pw
  # 
  # pdf(paste0(fig_path,"Figure1__3.pdf" ))
  # plot(pw)
  # dev.off()
 

c_vec_long =seq(0, 1.2, length.out = 40)

df_post_MTM_alpha(c_vec = c_vec_long, post= final,  alpha_l, N = 215)




################################
######## Univariate DMP ########
################################

data_raw = loadRData("~/Documents/GitHub/BNPconsistency/data/Slc.RData")
data<- list()
data$y = unname(data_raw)
alpha = 1
K_ = 10
M_it = 50000
nburn = 30000
N = dim(data_raw)[1]

alpha_l = c(0.01, 0.5,1, 2)
# final = comparison_data_alpha(data= data,K_= 10 , M_it = M_it, nburn= nburn, alpha_l=alpha_l)
# save(final, file = "~/Documents/GitHub/BNPconsistency/data/df_Slc.RData")

final = loadRData( "~/Documents/GitHub/BNPconsistency/data/df_Slc.RData")

## Normal plots

# fig_path = "../figures/Figure1/"
fig_df <- final

fig_df_mut <- fig_df$line%>%group_by(Process_type,Al)%>%mutate(pkn =density/sum(density))


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
                  Pkn_4 = round(julia_eval("Pkn_Dirichlet_mult.(1:10,215, 10, A4)"),3))%>% gather(Process_type, pkn,Pkn_1:Pkn_4)

df_prior$Type =rep("Prior", dim(df_prior)[1])
fig_df_mut$Type= rep("Posterior", dim(fig_df_mut)[1])

df_merged = rbind(df_prior,fig_df_mut[, c("K","Process_type", "pkn", "Type")] )


p_h<- ggplot(df_merged, aes(K,pkn,color =Process_type))+geom_bar(aes(linetype=Type),size = 0.7, stat="identity",alpha =0.0, position = "identity")+
  geom_vline(xintercept=2,  linetype="dashed")+ylab("Density")+xlab(TeX('$K_n$'))+
  theme_minimal()+theme(strip.background = element_blank(), strip.text.x = element_blank())+
  scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D",name = TeX(sprintf('$\\bar{\\alpha}$')), labels=unname(TeX(c(sprintf('$\\bar{\\alpha}$ = %.2f', alpha_l[1]),
                                                                                                                                                   sprintf('$\\bar{\\alpha}$ = %.1f', alpha_l[2]),
                                                                                                                                                   sprintf('$\\bar{\\alpha}$ = %1.f', alpha_l[3]),
                                                                                                                                                   sprintf('$\\bar{\\alpha}$ = %1.f', alpha_l[4])))))+
  scale_x_continuous(breaks = c(1,2,4,6,8,10), limits = c(0,11))+
  scale_linetype_manual(name = "Distribution",values = c(1, 2),guide = guide_legend(override.aes = list(linetype = c(1, 2),color = "black") ) )+
  facet_wrap(~Process_type)
p_h
pdf(paste("Figure1_Slc_DMP.pdf" ))
plot(p_h)
dev.off()


c_list = list(seq(0, 1e-6, length.out = 40), seq(0, 0.1, length.out = 40), seq(0, 0.12, length.out = 40), seq(0, 0.11, length.out = 40))

df_mut = plt_univ_MTM(c_list = c_list, post= final,  alpha_l, N = N)

pm <- ggplot(df_mut, aes(x = c, y = val, color = alpha))+ geom_line(aes(linetype=type))+
  scale_linetype_manual(name = "Distribution",values = c("solid", "dotdash"),guide = guide_legend(override.aes = list(color = "black")))+
  scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D",name = TeX(sprintf('$\\bar{\\alpha}$')), labels=unname(TeX(c(sprintf('$\\bar{\\alpha}$ = %.2f', alpha_l[1]),
                                                                                                                                       sprintf('$\\bar{\\alpha}$ = %.1f', alpha_l[2]),
                                                                                                                                       sprintf('$\\bar{\\alpha}$ = %.1f', alpha_l[3]),
                                                                                                                                       sprintf('$\\bar{\\alpha}$ = %1.f', alpha_l[4])))))+
  theme_minimal()+theme(strip.background = element_blank(), strip.text.x = element_blank())+
  ylab(expression(tilde(K)))+xlab(TeX('$c$')) +geom_hline(yintercept = 2, linetype="dashed", size = 0.5, alpha =0.5) +    facet_wrap(~alpha, scales="free_x")
plot(pm)
pdf(paste0("Figure_SLC_DMP_MTM.pdf"))
plot(pm)
dev.off()