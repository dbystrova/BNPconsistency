
require(tidyr)
require(e1071)
require(MASS)
require(MCMCpack)
require(mvtnorm)
require(Runuran)
require(flexclust)
library(cowplot)
library(ggplot2)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

l2norm_MV_gaussian = function(theta.j, theta.i){
  l2_mu = sum((theta.j[[1]] -theta.i[[1]] )^2)
  #l2_Sigma = sqrt(sum((theta.j[[2]] -theta.i[[2]] )^2))
  #l2_Sigma =norm(theta.j[[2]] -theta.i[[2]],"2")^2
  #l2_Sigma = sum(((theta.j[[2]] - theta.i[[2]])[upper.tri(theta.j[[2]] -theta.i[[2]], diag = TRUE)])^2)
  l2_Sigma =norm(theta.j[[2]] - theta.i[[2]],"F")^2
  return(sqrt(l2_mu + l2_Sigma ))
}

MTM<- function(G, w_n, c){
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
        if(l2norm_MV_gaussian(theta.k[[tau[j]]], theta.k[[tau[i]]]) <= w_n){
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
  
th =theta.k[tau]
p = p_new[tau]

sorted_weights  = sort(p, decreasing = TRUE, index.return=TRUE)
p_sorted = p[sorted_weights$ix]
th_sorted = th[sorted_weights$ix]
#stage 2



A = which(p_sorted > (c*w_n)^2)
N =  which(p_sorted <= (c*w_n)^2)

for(i in A){
  for(j in A){
    if (j < i){
    #  print(p[i]*(l2norm_MV_gaussian(th[[i]],th[[j]]))^2)
      if (p_sorted[i]*(l2norm_MV_gaussian(th_sorted[[i]],th_sorted[[j]]))^2 <= (c*w_n)^2){
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
    if (v<= l2norm_MV_gaussian(th_sorted[[j]],th_sorted[[i]])){
    v = l2norm_MV_gaussian(th_sorted[[j]],th_sorted[[i]])
    ind_min = j
   }
  }
  p_sorted[ind_min] = p_sorted[ind_min] + p_sorted[i] 
 }
return (list(p=p_sorted[A], th=th_sorted[A]))
}


apply_MTM<- function(p_eta, p_mu, p_sig,c,n, w_n = (log(n)/n)^(1/4)){
  post_k<-c()
  it <- dim(p_eta)[1]
  for (k in 1:it ){
    G <- list()
    G$p = p_eta[k,]
    G$theta = list ()
    for (i in 1:length(p_eta[k,])){
      G$theta[[i]]<-  list(p_mu[k,,i],p_sig[k,,,i])
    }
    G.post <- MTM(G, w_n, c)
    #print(G.post$p)
    post_k[k]<- length(G.post$p)
  }
  
  return(post_k)
}


df_post<- function(n, c_vec,ind, post){
  df_post_k= tibble(it= 1:dim(post$eta[[1]])[1],
                    Pkn_c1 = round(apply_MTM(post$eta[[ind]],post$mu[[ind]],post$sig[[ind]], c =c_vec[1], n=n),3),
                    Pkn_c2= round(apply_MTM(post$eta[[ind]],post$mu[[ind]],post$sig[[ind]], c = c_vec[2], n=n),3), 
                    Pkn_c3 = round(apply_MTM(post$eta[[ind]],post$mu[[ind]],post$sig[[ind]], c =c_vec[3], n=n),3),
                    Pkn_c4 = round(apply_MTM(post$eta[[ind]],post$mu[[ind]],post$sig[[ind]], c = c_vec[4], n=n),3))%>% gather(Process_type, pkn,Pkn_c1:Pkn_c4)
  
  df_post_k$N = rep(n,dim(df_post_k)[1])
  return(df_post_k)
  
}
