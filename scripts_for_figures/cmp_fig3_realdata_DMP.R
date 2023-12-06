source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Random_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Estimation_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Identification_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Gibbs_sampling_function.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/plt_Figure1.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")

l2norm_MV_gaussian = function(theta.j, theta.i, diag = FALSE){
  l2_mu <- sum((theta.j[[1]] -theta.i[[1]] )^2)
  if (diag){l2_Sigma <- sum((theta.j[[2]] -theta.i[[2]] )^2) }
  else{
    l2_Sigma <- norm(theta.j[[2]] - theta.i[[2]],"F")^2}
  return(sqrt(l2_mu + l2_Sigma))
}

l2norm_univ <- function(theta.j, theta.i){
  l2_mu <- sum((theta.j[[1]] -theta.i[[1]])^2)
  l2_Sigma <- sum((theta.j[[2]] -theta.i[[2]])^2)
  return(sqrt(l2_mu + l2_Sigma))
}

apply_MTM_univ <- function(p_eta, p_mu, p_sig, c, n, w_n = (log(n)/n)^(1/4)){
  post_k<-c()
  it <- dim(p_eta)[1]
  for (k in 1:it){
    G <- list()
    G$p <- p_eta[k,]
    G$theta <- list ()
    for (i in 1:length(p_eta[k,])){
      G$theta[[i]] <-  list(p_mu[k,i], p_sig[k,i])
    }
    G.post <- MTM_univ(G, w_n, c)
    #print(G.post$p)
    post_k[k] <- length(G.post$p)
  }
  
  return(post_k)
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

df_post_MTM_alpha<- function(c_vec, post, alpha_l, N, name){
  n_it <- dim( post$eta[[1]])[1]
  tokeep <- seq(2, n_it, 4)
  for  (i in 1:length(alpha_l)){
    post$eta[[i]] <- post$eta[[i]][tokeep,]
    post$mu[[i]] <- post$mu[[i]][tokeep,,] 
    post$sig[[i]] <- post$sig[[i]][tokeep,,,] 
  }
  df_post_k <- tibble(с=c_vec)
  df_post_k_map <- tibble(с=c_vec)
  for (i in 1:length(c_vec)){
    for (j in 1:length(alpha_l)){
      result_MTM <- apply_MTM(post$eta[[j]], post$mu[[j]], post$sig[[j]], c=c_vec[i], n=N)
      df_post_k[i,(j+1)] <- mean(round(result_MTM, 3))
      map <- names(table(result_MTM)[which.max(table(result_MTM))])
      df_post_k_map[i,(j+1)] <- as.numeric(map)
    }
  }
  names(df_post_k) <- c('c', as.character(alpha_l))
  df_mut_mean <- df_post_k %>% gather(key=alpha_level, value=val, 2:5)
  df_mut_mean$type <- rep("Mean", dim(df_mut_mean)[1])
  names(df_post_k_map) <- c('c', as.character(alpha_l))
  df_mut_map <- df_post_k_map %>% gather(key=alpha_level, value=val, 2:5)
  df_mut_map$type <- rep("MAP", dim(df_mut_map)[1])
  df_mut <- rbind(df_mut_mean,df_mut_map)
  df_mut$alpha_level <- as.factor(df_mut$alpha_level)
  
  return(df_mut)
}

df_post_univ <- function(c_vec, post, ind, alpha, N){
  df_post_k <- tibble(с_val=c_vec)
  df_post_k$val <- rep(NA, length(c_vec))
  df_post_k_map <- tibble(с_val=c_vec)
  df_post_k_map$val <- rep(NA, length(c_vec))
  
  n_it <- dim( post$eta[[1]])[1]
  tokeep <- seq(2, n_it, 4)
  post$eta[[ind]] <- post$eta[[ind]][tokeep,]
  post$mu[[ind]] <- post$mu[[ind]][tokeep,,] 
  post$sig[[ind]] <- post$sig[[ind]][tokeep,,,] 
  
  for (i in 1:length(c_vec)){
    result_MTM <- apply_MTM_univ(post$eta[[ind]], post$mu[[ind]], post$sig[[ind]], c=c_vec[i], n=N)
    map <- names(table(result_MTM)[which.max(table(result_MTM))])
    df_post_k_map$val[i] <- as.numeric(map)
    df_post_k$val[i] <- mean(round(result_MTM,3))
  }
  #df_mut = df_post_k%>%gather(Process_type, pkn,2:21)
  df_post_k$alpha <- rep(alpha, length(c_vec))
  df_post_k$type <- rep("Mean", length(c_vec))
  df_post_k_map$alpha <- rep(alpha, length(c_vec))
  df_post_k_map$type <- rep("MAP", length(c_vec))
  return(list(df_post_k=df_post_k, df_post_k_map=df_post_k_map))
}

cmp_univ_MTM<- function(c_list, post,  alpha_l, N){
  df_1 <- df_post_univ(c_list[[1]], post, ind=1, alpha=alpha_l[1], N)
  df_2 <- df_post_univ(c_list[[2]], post, ind=2, alpha=alpha_l[2], N)
  df_3 <- df_post_univ(c_list[[3]], post, ind=3, alpha=alpha_l[3], N)
  df_4 <- df_post_univ(c_list[[4]], post, ind=4, alpha=alpha_l[4], N)
  
  df_fin <- rbind(df_1$df_post_k, df_2$df_post_k, df_3$df_post_k, df_4$df_post_k)
  df_fin_map <- rbind(df_1$df_post_k_map, df_2$df_post_k_map, df_3$df_post_k_map, df_4$df_post_k_map)
  
  df_mut <- rbind(df_fin_map, df_fin)
  df_mut$alpha <- as.factor(df_mut$alpha)
  
  return(df_mut)
}

##################################
######## Multivariate DMP ########
##################################

final_mult <- loadRData( "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_thyroid_DMP.RData")
alpha_l <- as.numeric(levels(as.factor(final_mult$line$Al)))
N <- final_mult$hist$N[1]

c_vec_long <- seq(0, 1.2, length.out = 40)

df <- df_post_MTM_alpha(c_vec = c_vec_long, post = final_mult, alpha_l, N = N)
save(df, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_thyroid_DMP_MTM.RData")

################################
######## Univariate DMP ########
################################

final_univ <- loadRData( "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_Slc_DMP.RData")
alpha_l <- as.numeric(levels(as.factor(final_univ$line$Al)))
N <- final_univ$hist$N[1]

c_list <- list(seq(0, 1e-6, length.out = 40), seq(0, 0.1, length.out = 40), 
               seq(0, 0.12, length.out = 40), seq(0, 0.11, length.out = 40))

df <- cmp_univ_MTM(c_list = c_list, post = final_univ, alpha_l, N = N)
save(df, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_Slc_DMP_MTM.RData")
