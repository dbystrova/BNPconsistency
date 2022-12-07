## read sources
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Random_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Estimation_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Identification_SpMix.R")

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
library(abind)
#---------- B) Specification of the simulation and prior parameters -----------------------------------------------

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

compute_log_lik<- function(K, y, M, burnin, Eta, Mu, Sigma){
  n_data = dim(y)[1]
  plk = 0
  lik = matrix(0,n_data,K )
  log_lik = c()
  for (i in 1:(M-burnin)){
    for (j in 1:K){
      lik[,j] = Eta[burnin+i,j]*dmvnorm(y,Mu[burnin+i,,j], Sigma[burnin+i,,,j])
      #k_ = S_alt_matrix[j,i]
      #p_l <- sum(sapply(1:n_data, function(i) Eta[j,S_alt_matrix[j,i]] * dmvnorm(y[i,],Mu[j,,S_alt_matrix[j,i]], Sigma[j,,,S_alt_matrix[j,i]])))
      #k_ = S_alt_matrix[j,i]
      #plk = plk + Eta[j,k_] * dmvnorm(y[i,], Mu[j,,k_], Sigma[j,,,k_])
    }
    log_lik[i] = sum(log(apply(lik,1,sum)))
  }
  return(log_lik)
  }



MCMC_function <- function(data, e0=0.01, K, M, burnin,seed_ = 1000, priorOnE0 = FALSE) {
  print(seed_)
  y <- as.matrix(data$y)
  Mmax <- M + burnin
  
  ## read dimensions of data:
  r <- length(y[1, ])
  N <- length(y[, 1])
  
  ## Dirichlet parameter for the mixture weights
  ## variance of the normal proposal for the MH step for estimating e0
  c_proposal <- 0.8
  
  R <- apply(y, 2, function(x) diff(range(x)))
  
  ## prior on Sigma_k
  c0 <- 2.5 + (r - 1)/2
  C0 <- 0.75 * cov(y)
  g0 <- 0.5 + (r - 1)/2
  G0 <- 100 * g0/c0 * diag((1/R^2))
  
  ## prior on mu
  b0 <- apply(y, 2, median)
  B_0 <- rep(1, r)  #initial values for lambda are 1
  B0 <- diag((R^2) * B_0)
  nu <- 0.5
  
  ## initial values for parameters to be estimated:
  eta_0 <- rep(1/K, K)
  sigma_0 <- array(0, dim = c(r, r, K))
  for (k in 1:K) {
    sigma_0[, , k] <- C0
  }
  C0_0 <- C0
  
  
  ## initial classification
  groups <- K
  cl_y <- kmeans(y, centers = groups, nstart = 30)
  S_0 <- cl_y$cluster
  mu_0 <- cbind(t(cl_y$centers))
  
  
  ## generate matrices for saving the results:
  Eta <- matrix(0, M, K)
  Mu <- array(0, dim = c(M, r, K))
  B <- matrix(0, M, r)
  Eta_Matrix_FS <- matrix(0, r, K)
  Sigma_Matrix_FS <- array(0, dim = c(r, r, K))
  Mu_Matrix_FS <- array(0, dim = c(r, K))
  
  
  
  #---------- C) Gibbs sampling from the posterior -----------------------------------------------
  #print(B_0)
  ################ call MCMC procedure
  ############### run 2 parallel chains 
  estGibbs <- MultVar_NormMixt_Gibbs_IndPriorNormalgamma(y, S_0, mu_0, sigma_0, eta_0, e0, c0, C0_0, 
                                                         g0, G0, b0, B0, nu, B_0, M, burnin, c_proposal, priorOnE0 = priorOnE0, lambda = FALSE,seed =seed_)
  estGibbs_2 <- MultVar_NormMixt_Gibbs_IndPriorNormalgamma(y, S_0, mu_0, sigma_0, eta_0, e0, c0, C0_0, 
                                                         g0, G0, b0, B0, nu, B_0, M, burnin, c_proposal, priorOnE0 = priorOnE0, lambda = FALSE,seed =seed_+1)
  
  
  Mu <-estGibbs$Mu
  Sigma <- estGibbs$Sigma
  Eta <- estGibbs$Eta
  S_alt_matrix <- estGibbs$S_alt_matrix
  B <- estGibbs$B
  e0_vector <- estGibbs$e0_vector
  acc_rate <- estGibbs$acc_rate
        
  Eta_combined =   rbind(estGibbs$Eta[(burnin+1):M,],estGibbs_2$Eta[(burnin+1):M,])   
  Mu_combined =   abind(estGibbs$Mu[(burnin+1):M,,],estGibbs_2$Mu[(burnin+1):M,,] ,along = 1 )                                
  Sigma_combined =   abind(estGibbs$Sigma[(burnin+1):M,,,],estGibbs_2$Sigma[(burnin+1):M,,,] ,along = 1 )                                
  
  Nk_matrix_alt <- estGibbs$Nk_matrix_alt
  Nk_matrix_alt2 <- estGibbs_2$Nk_matrix_alt
  nonnormpost_mode_list <- estGibbs$nonnormpost_mode_list

  ## comute log-likelihood for each chain
  ll1<- compute_log_lik(K, y, M, burnin, Eta, Mu, Sigma)
  ll2<- compute_log_lik(K, y, M, burnin, estGibbs_2$Eta,estGibbs_2$Mu,estGibbs_2$Sigma)
  ## covergence diagnostics
  log_lik_combines <- mcmc.list(mcmc(ll1),mcmc(ll2))
  Rhat_ll<- gelman.diag(log_lik_combines)$psrf[1]
  
  ##### number of nonempty components (nne_gr), after burnin:
  #nne_gr <- apply(Nk_matrix_alt, 1, function(x) sum(x != 0))
  #table(nne_gr)
  #plot(1:length(nne_gr), nne_gr, type = "l", ylim = c(1, max(nne_gr)), xlab = paste("iteration"), main = "number of non-empty components")
  
  
  
  
  
  ## convergence diagnostic
  #Sigma_es = apply(Sigma, c(2,3,4),effectiveSize)
  #Sigma_es_mean = mean(Sigma_es)
  #Mu_es = apply(Mu, c(2,3),effectiveSize)
  #Mu_es_mean = mean(Mu_es)
  ##### number of nonempty components (nne_gr), after burnin:
  nne_gr <- apply(Nk_matrix_alt2, 1, function(x) sum(x != 0))
  table(nne_gr)
  #plot(1:length(nne_gr), nne_gr, type = "l", ylim = c(1, max(nne_gr)), 
  #     xlab = paste("iteration"), main = "number of non-empty components")
  
  #---------- E) Identification of the mixture model -----------------------------------------------
  ##### estimating the number of nonempty components
  K0_vector <- rowSums(Nk_matrix_alt[(burnin+1):M,] != 0)  #vector with number of non-empty groups of each iteration
  K0_vector2 <- rowSums(Nk_matrix_alt2[(burnin+1):M,] != 0)  #vector with number of non-empty groups of each iteration
  p_K0 <- tabulate(c(K0_vector,K0_vector2), K)
  p_K0
  #par(mfrow = c(1, 1))
  #barplot(p_K0, names = 1:K, xlab = "number of non-empty groups K0", col = "green", ylab = "freq")
  K0 <- which.max(p_K0)
  K0  #mode K0 is the estimator for K_true
  M0 <- sum(K0_vector == K0)
  M0  #M0 draws have exactly K0 non-empty groups
  
return(list(p_k = p_K0,ll_rhat= Rhat_ll,p_k_all = c(K0_vector,K0_vector2), Eta = Eta_combined, Mu= Mu_combined, Sigma = Sigma_combined ))
}





MCMC_function_ra <- function(data, e0=0.01, K, M, burnin,seed_ = 1000) {
  print(seed_)
  y <- as.matrix(data$y)
  Mmax <- M + burnin
  
  ## read dimensions of data:
  r <- length(y[1, ])
  N <- length(y[, 1])
  
  ## Dirichlet parameter for the mixture weights
  ## variance of the normal proposal for the MH step for estimating e0
  c_proposal <- 0.8
  
  R <- apply(y, 2, function(x) diff(range(x)))
  
  ## prior on Sigma_k
  c0 <- 2.5 + (r - 1)/2
  C0 <- 0.75 * cov(y)
  g0 <- 0.5 + (r - 1)/2
  G0 <- 100 * g0/c0 * diag((1/R^2))
  
  ## prior on mu
  b0 <- apply(y, 2, median)
  B_0 <- rep(1, r)  #initial values for lambda are 1
  B0 <- diag((R^2) * B_0)
  nu <- 0.5
  
  ## initial values for parameters to be estimated:
  eta_0 <- rep(1/K, K)
  sigma_0 <- array(0, dim = c(r, r, K))
  for (k in 1:K) {
    sigma_0[, , k] <- C0
  }
  C0_0 <- C0
  
  
  ## initial classification
  groups <- K
  cl_y <- kmeans(y, centers = groups, nstart = 30)
  S_0 <- cl_y$cluster
  mu_0 <- cbind(t(cl_y$centers))
  
  
  ## generate matrices for saving the results:
  Eta <- matrix(0, M, K)
  Mu <- array(0, dim = c(M, r, K))
  B <- matrix(0, M, r)
  Eta_Matrix_FS <- matrix(0, r, K)
  Sigma_Matrix_FS <- array(0, dim = c(r, r, K))
  Mu_Matrix_FS <- array(0, dim = c(r, K))
  
  
  
  #---------- C) Gibbs sampling from the posterior -----------------------------------------------
  #print(B_0)
  ################ call MCMC procedure
  ############### run 2 parallel chains 
  estGibbs <- MultVar_NormMixt_Gibbs_IndPriorNormalgamma(y, S_0, mu_0, sigma_0, eta_0, e0=1, c0, C0_0, 
                                                         g0, G0, b0, B0, nu, B_0, M, burnin, c_proposal, priorOnE0 = TRUE, lambda = FALSE,seed =seed_)
  estGibbs_2 <- MultVar_NormMixt_Gibbs_IndPriorNormalgamma(y, S_0, mu_0, sigma_0, eta_0, e0=1, c0, C0_0, 
                                                           g0, G0, b0, B0, nu, B_0, M, burnin, c_proposal, priorOnE0 = TRUE, lambda = FALSE,seed =seed_+1)
  
  
  Mu <-estGibbs$Mu
  Sigma <- estGibbs$Sigma
  Eta <- estGibbs$Eta
  S_alt_matrix <- estGibbs$S_alt_matrix
  B <- estGibbs$B
  e0_vector <- estGibbs$e0_vector
  acc_rate <- estGibbs$acc_rate
  
  Eta_combined =   rbind(estGibbs$Eta[(burnin+1):M,],estGibbs_2$Eta[(burnin+1):M,])   
  Mu_combined =   abind(estGibbs$Mu[(burnin+1):M,,],estGibbs_2$Mu[(burnin+1):M,,] ,along = 1 )                                
  Sigma_combined =   abind(estGibbs$Sigma[(burnin+1):M,,,],estGibbs_2$Sigma[(burnin+1):M,,,] ,along = 1 )                                
  e0_combined =   c(estGibbs$e0_vector[(burnin+1):M],estGibbs_2$e0_vector[(burnin+1):M])                                
  
  Nk_matrix_alt <- estGibbs$Nk_matrix_alt
  Nk_matrix_alt2 <- estGibbs_2$Nk_matrix_alt
  nonnormpost_mode_list <- estGibbs$nonnormpost_mode_list
  
  ## comute log-likelihood for each chain
  ll1<- compute_log_lik(K, y, M, burnin, Eta, Mu, Sigma)
  ll2<- compute_log_lik(K, y, M, burnin, estGibbs_2$Eta,estGibbs_2$Mu,estGibbs_2$Sigma)
  ## covergence diagnostics
  log_lik_combines <- mcmc.list(mcmc(ll1),mcmc(ll2))
  Rhat_ll<- gelman.diag(log_lik_combines)$psrf[1]
  ##### number of nonempty components (nne_gr), after burnin:
  nne_gr <- apply(Nk_matrix_alt2, 1, function(x) sum(x != 0))
  table(nne_gr)
  #plot(1:length(nne_gr), nne_gr, type = "l", ylim = c(1, max(nne_gr)), 
  #     xlab = paste("iteration"), main = "number of non-empty components")
  
  #---------- E) Identification of the mixture model -----------------------------------------------
  ##### estimating the number of nonempty components
  K0_vector <- rowSums(Nk_matrix_alt[(burnin+1):M,] != 0)  #vector with number of non-empty groups of each iteration
  K0_vector2 <- rowSums(Nk_matrix_alt2[(burnin+1):M,] != 0)  #vector with number of non-empty groups of each iteration
  p_K0 <- tabulate(c(K0_vector,K0_vector2), K)
  p_K0
  #par(mfrow = c(1, 1))
  #barplot(p_K0, names = 1:K, xlab = "number of non-empty groups K0", col = "green", ylab = "freq")
  K0 <- which.max(p_K0)
  K0  #mode K0 is the estimator for K_true
  M0 <- sum(K0_vector == K0)
  M0  #M0 draws have exactly K0 non-empty groups
  
  return(list(p_k = p_K0,ll_rhat= Rhat_ll,p_k_all = c(K0_vector,K0_vector2), Eta = Eta_combined, Mu= Mu_combined, Sigma = Sigma_combined, e0 =e0_combined ))
}


comparison_n<- function(ds_list, alpha,K_, M_it, nburn){
  pk<- list()
  N<- c()
  R_h <- c()
  W_non_sorted <- list()
  W <- list()
  Mu_mat <- list()
  S_mat<- list()
  for (i in 1:length(ds_list)){
    data =  loadRData(ds_list[i])
    pk[[i]] <- MCMC_function(data, e0=alpha, K=K_, M=M_it, burnin=nburn)
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
  
  df$alpha = c(rep(alpha,dim(df_)[1]*length(ds_list))) 
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
  df_l$alpha = c(rep(alpha,(M_it -nburn)*2*length(ds_list))) 
  df_l$Rh = df_l2$Rh_val
  df_l$N =df_l3$N_val
  return(list(line = df, hist =df_l, weights = df4_post_, eta = W_non_sorted, mu = Mu_mat, sig = S_mat))
}



comparison_n_2<- function(ds_list,K_, M_it, nburn, alpha_l){
  pk<- list()
  N<- c()
  R_h <- c()
  W_non_sorted <- list()
  W <- list()
  Mu_mat <- list()
  S_mat<- list()
  for (i in 1:length(ds_list)){
    data =  loadRData(ds_list[i])
    pk[[i]] <- MCMC_function(data, e0=alpha_l[i]/K_, K=K_, M=M_it, burnin=nburn)
    N[i]<- dim(data$y)[1]
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
  df5_ = tibble(K= 1:K_)
  for (j in 1:length(ds_list)){
    name_ <- paste("Pkn_", j, sep = "")
    name2_ <- paste("Rh_", j, sep = "")
    name3_ <- paste("N_", j, sep = "")
    name4_ <- paste("Alpha_", j, sep = "")
    name5_ <- paste("W_", j, sep = "")
    df_[,name_]<- pk[[j]]$p_k
    df2_[,name2_]<- rep(pk[[j]]$ll_rhat,length(pk[[j]]$p_k))
    df3_[,name3_]<- rep(N[j],length(pk[[j]]$p_k))
    df4_[,name4_]<- rep(alpha_l[j]/K_,length(pk[[j]]$p_k))
    df5_[,name5_]<- rep(N[j],length(pk[[j]]$p_k))
  }
  df = df_%>% gather(Process_type, density,  paste("Pkn_", 1, sep = ""):paste("Pkn_", length(ds_list), sep = ""))
  df2 = df2_%>% gather(Rh, Rh_val,  paste("Rh_", 1, sep = ""):paste("Rh_", length(ds_list), sep = ""))
  df3 = df3_%>% gather(N_, N_val,  paste("N_", 1, sep = ""):paste("N_", length(ds_list), sep = ""))
  df4 = df4_%>% gather(Alpha_, Alpha_val,  paste("Alpha_", 1, sep = ""):paste("Alpha_", length(ds_list), sep = ""))
  df5 = df5_%>% gather(W_, W_val,  paste("W_", 1, sep = ""):paste("W_", length(ds_list), sep = ""))
  
  W_df <- do.call(cbind, W)
  df5_post <- cbind(df5,t(W_df))
  
  df5_post_ <- gather(df5_post, key = "it",value="weights", 4: dim(df5_post)[2])
  
  
  
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
  df_l$Rh = df_l2$Rh_val
  df_l$N =df_l3$N_val
  df_l$Al =df_l4$Alpha_val
  return(list(line = df, hist =df_l, weights = df5_post_, eta = W_non_sorted, mu = Mu_mat, sig = S_mat))
}

