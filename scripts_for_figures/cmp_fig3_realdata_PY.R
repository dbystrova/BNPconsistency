library(BNPmix)
library(tidyverse)
library(latex2exp)
library(viridis)
library(MCMCpack)
library(JuliaCall)
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")

l2norm_univ <- function(theta.j, theta.i){
  l2_mu <- sum((theta.j[[1]] -theta.i[[1]])^2)
  l2_Sigma <- sum((theta.j[[2]] -theta.i[[2]])^2)
  return(sqrt(l2_mu + l2_Sigma))
}

MTM_univ <- function(G, w_n, c){
  p.k <- G$p
  theta.k <- G$theta
  ## stage 1
  if (length(which(p.k == 0)) > 0) {
    theta.k <- theta.k[-c(which(p.k == 0))]
    p.k <- p.k[-c(which(p.k == 0))]
  }
  theta_ind <- sample(1:length(theta.k), length(theta.k), replace=FALSE, prob=p.k)
  tau <- theta_ind
  
  p_new <- vector(mode="numeric", length=length(p.k ))
  p_new <- p.k
  t <- length(tau)
  i <- 1 
  while (i<=t){
    j <- 1
    while (j<=t){
      if(tau[j]<tau[i]){
        if(l2norm_univ(theta.k[[tau[j]]], theta.k[[tau[i]]]) <= w_n){
          p_new[tau[i]] <- p.k[tau[j]] + p.k[tau[i]]
          tau <- tau[-j]
          t <- length(tau)
          if (j<i) {i <- i-1}
          j <- j-1
        }
      }
      j <- j + 1
    }
    i <- i +1
  }
  
  #stage 2
  sorted_weights <- sort(p_new[tau], decreasing = TRUE, index.return=TRUE)
  p <- sorted_weights$x
  th <- theta.k[tau]
  th_sorted <- th[sorted_weights$ix]
  
  A <- which(p > (c*w_n)^2)
  N <- which(p <= (c*w_n)^2)
  
  for(i in A){
    for(j in A){
      if (j < i){
        if (p[i]*(l2norm_univ(th[[i]],th[[j]]))^2 <= (c*w_n)^2){
          N <- c(N, i)
          A <- A[-i]
        }
      }
    }
  }
  
  for(i in N){
    v <- 0
    ind_min <- 1
    for (j in A){
      if (v <= l2norm_univ(th[[j]],th[[i]])){
        v <- l2norm_univ(th[[j]],th[[i]])
        ind_min <- j
      }
    }
    p[ind_min] <- p[ind_min] + p[i] 
  }
  return (list(p=p[A], th=th[A]))
}


apply_MTM_BNPmix<- function(fit_obj, c, n,  w_n = (log(n)/n)^(1/4), it = 10000){
  post_k <- list()
  it <- dim(fit_obj$clust)[1]
  n <- dim(fit_obj$clust)[2]
  for (k in seq(2, it, 4)){
    G <- list()
    G$p <- fit_obj$probs[[k]]
    G$theta <- list()
    for (i in 1:length(fit_obj$probs[[k]])){
      G$theta[[i]] <-  list(fit_obj$mean[[k]][i,], fit_obj$sigma2[[k]][,,i]) 
    }
    if(length(G$p) == length(G$theta))
      G.post <- MTM(G, w_n, c)
    # else
    #   cat("iteration :", k, ", weights : ", length(G$p), ", atoms : ", length(G$theta), "\n")
    post_k <- append(post_k, length(G.post$p))
  }
  
  return(post_k)
}

apply_MTM_BNPmix_univ <- function(fit_obj, c, n,  w_n = (log(n)/n)^(1/4)){
  post_k <- list()
  it <- dim(fit_obj$clust)[1]
  n <- dim(fit_obj$clust)[2]
  for (k in seq(2, it, 4)){
    G <- list()
    G$p <- c(fit_obj$probs[[k]])
    G$theta <- list()
    for (i in 1:length(fit_obj$mean[[k]])){
      G$theta[[i]] <-  list(fit_obj$mean[[k]][i,], fit_obj$sigma2[[k]][i,]) 
    }
    if(length(G$p) == length(G$theta))
      G.post <- MTM_univ(G, w_n, c)
    # else
    #   cat("iteration :", k, ", weights : ", length(G$p), ", atoms : ", length(G$theta), "\n")
    post_k <- append(post_k, length(G.post$p))
  }
  
  return(post_k)
}

py_mtm <- function(final, c_vec){
  fit.sim1 <- final$f1
  fit.sim2 <- final$f2
  fit.sim3 <- final$f3
  fit.sim4 <- final$f4
  
  mean <- rep(0,length(c_vec))
  map <- rep(0,length(c_vec))
  n <- final$N
  if(final$univ){
    for(i in 1:length(c_vec)){
      result_MTM1 <- unlist(apply_MTM_BNPmix_univ(fit.sim1, c = c_vec[i], w_n = (log(n))^(-2)))
      result_MTM2 <- unlist(apply_MTM_BNPmix_univ(fit.sim2, c = c_vec[i], w_n = (log(n))^(-2)))
      result_MTM3 <- unlist(apply_MTM_BNPmix_univ(fit.sim3, c = c_vec[i], w_n = (log(n))^(-2)))
      result_MTM4 <- unlist(apply_MTM_BNPmix_univ(fit.sim4, c = c_vec[i], w_n = (log(n))^(-2)))
      result_MTM <- c(result_MTM1, result_MTM2, result_MTM3, result_MTM4)
      map[i] <- as.numeric(names(table(result_MTM)[which.max(table(result_MTM))]))
      mean[i] <- mean(round(result_MTM))
    }
  }
  else{
    for(i in 1:length(c_vec)){
      result_MTM1 <- unlist(apply_MTM_BNPmix(fit.sim1, c = c_vec[i], w_n = (log(n))^(-2)))
      result_MTM2 <- unlist(apply_MTM_BNPmix(fit.sim2, c = c_vec[i], w_n = (log(n))^(-2)))
      result_MTM3 <- unlist(apply_MTM_BNPmix(fit.sim3, c = c_vec[i], w_n = (log(n))^(-2)))
      result_MTM4 <- unlist(apply_MTM_BNPmix(fit.sim4, c = c_vec[i], w_n = (log(n))^(-2)))
      result_MTM <- c(result_MTM1, result_MTM2, result_MTM3, result_MTM4)
      map[i] <- as.numeric(names(table(result_MTM)[which.max(table(result_MTM))]))
      mean[i] <- mean(round(result_MTM))
    }
  }
  MTM_out <- tibble(c = c_vec, mean = mean, map = map, 
                   alpha = rep(final$Al, length(c_vec)), sigma = rep(final$Sig, length(c_vec)))
  return(MTM_out)
}

#################################
######## Multivariate PY ########
#################################
c_vec = seq(0, 13, length.out = 50)

final <- loadRData("~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_thyroid_PY.RData")
for(i in 1:length(final)){
  fin_i <- final[[i]]
  if(i == 1){
    MTM_out <- py_mtm(fin_i, c_vec = c_vec)
    MTM_out$Process_type <- paste0("Pkn_",i)
  }
  else{
    MTM_out_ <- py_mtm(fin_i, c_vec = c_vec)
    MTM_out_$Process_type <- paste0("Pkn_",i)
    MTM_out <- rbind(MTM_out, MTM_out_)
  }
}
save(MTM_out, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_thyroid_PY_MTM.RData")

###############################
######## Univariate PY ########
###############################
c_list = list(seq(0, 3, length.out = 50), seq(0, 5, length.out = 50), seq(0, 5.5, length.out = 50), seq(0, 6, length.out = 50))

final <- loadRData("~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_Slc_PY.RData")
for(i in 1:length(final)){
  fin_i <- final[[i]]
  if(i == 1){
    MTM_out <- py_mtm(fin_i, c_vec = c_vec)
    MTM_out$Process_type <- paste0("Pkn_",i)
  }
  else{
    MTM_out_ <- py_mtm(fin_i, c_vec = c_vec)
    MTM_out_$Process_type <- paste0("Pkn_",i)
    MTM_out <- rbind(MTM_out, MTM_out_)
  }
}
save(MTM_out, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_Slc_PY_MTM.RData")
