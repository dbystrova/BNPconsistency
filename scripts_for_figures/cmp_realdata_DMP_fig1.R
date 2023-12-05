source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Random_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Estimation_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Identification_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Gibbs_sampling_function.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/plt_Figure1.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")
library(JuliaCall)

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

alpha = 1
K_ = 10
M_it = 50000
nburn = 30000

##################################
######## Multivariate DMP ########
##################################

data_raw = loadRData("~/Documents/GitHub/BNPconsistency/scripts_for_figures/real_data/thyroid.RData")
data<- list()
data$y = data_raw[,2:6]
N = dim( data_raw[,2:6])[1]

alpha_l = c(0.01,0.5, 1, 10)
final = comparison_data_alpha(data= data,K_= 10 , M_it = M_it, nburn= nburn, alpha_l=alpha_l)
save(final, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_thyroid_DMP.RData")



################################
######## Univariate DMP ########
################################

data_raw = loadRData("~/Documents/GitHub/BNPconsistency/scripts_for_figures/real_data/Slc.RData")
data<- list()
data$y = unname(data_raw)
N = dim(data_raw)[1]

alpha_l = c(0.01, 0.5,1, 2)
final = comparison_data_alpha(data= data,K_= 10 , M_it = M_it, nburn= nburn, alpha_l=alpha_l)
save(final, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_Slc_DMP.RData")
