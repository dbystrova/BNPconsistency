
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
fig6 <- loadRData("~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig6.RData")

G <- list()
G$p.k = 
G$theta = c(0.1,0.1)
c=1
n = 50 
w_n = 1/n

MTM<- function(G, w_n, c){
  p.k = G$p
  theta.k= G$theta
  
  ## stage 1
  new_theta = sample(theta.k, length(theta.k), replace=FALSE, prob = p.k)
  tau <- vector(mode="numeric", length=length(theta.k))
  for (i in 1:length(theta.k)){
    tau_k = which(new_theta ==theta.k[i] )
    tau[i ]= tau_k
  }
  print(tau)
  p_new = vector(mode="numeric", length=length(nj_k))
  p_new= p.k
  t = length(tau)
  i = 1 
  j = 1
while (i<t){
  while (j<t){
    if(tau[j]< tau[i]){
      l2norm = sqrt((theta_k[tau[j]] -  theta_k[tau[i]])^2 +(sigma_k[tau[j]] -  sigma_k[tau[i]])^2)
      if(l2norm <= w_n){
        p_new[tau[i]] =p.k[tau[j]] + p.k[tau[i]]
        p_new = p_new[- tau[j]]
        tau = tau[tau!=tau[j]]
        t = t-1
        print(c(i,j))
      }
    }
    j = j + 1
  }
  i = i +1
}

#stage 2
dat<- data.frame(p = p_new, th =theta_k[tau], sig=sigma_k[tau]  )
dat[order(dat$p.Freq, decreasing = TRUE),]
c = 0.5
A = which(dat$p.Freq > (c*w_n)^2)
N =  which(dat$p.Freq <= (c*w_n)^2)

t = length(N)
i=1
j=1

while ( i < t){
  while ( j < t){
    if (A[j] < A[i]){
      lnorm =(dat$th[A[i]]-dat$th[A[j]])^2 +(dat$sig[A[i]]-dat$sig[A[j]])^2
      if (dat$p.Freq[A[i]]*lnorm <= (c*w_n)^2){
        N= c(N, A[i])
        A=A[-i]
        t = t-1
      }
    }
  } 
  return ()
}


