
library(JuliaCall)
library(rust)
library(RcppArmadillo)
library(arm)
library(Rcpp)
library(gtools)
library(JuliaCall)

Rcpp::sourceCpp('~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/src/user_fns.cpp')
#based on the blogpost "https://www.kent.ac.uk/smsas/personal/msr/rlaptrans.html"
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/rlaptrans.r")



#test for ratio of uniform sampling 
julia <- julia_setup()
julia_library("GibbsTypePriors")
julia_library("DataFrames")
julia_library("DataFramesMeta")



lt.temp_st_pdf <- function(s, c, sigma, k) {
  exp( - c*( (s+k)^(sigma)  - k^(sigma) ))
}

pdf_lk_mat<- function(l,v, n_k, sigma,H, mat){
  return( (v^l)*exp(mat[n_k,l]))
}

sample_lk_mat<- function(nk_vec,v,sigma,H,M){
  l_post<-c()
  k<- length(nk_vec)
  for (i in 1:k){
    l_vec<- 1:nk_vec[i]
    if (length(l_vec)==1){
      l_post[i]=l_vec
    }
    else{
      p_v<- sapply(l_vec, function(x) pdf_lk_mat(x,v,nk_vec[i],sigma,H,mat=M))
      pv_norm<- p_v/sum(p_v)
      l_post[i]<- sample(1:(nk_vec[i]),size=1, replace=TRUE, prob=pv_norm)
    }
  }
  return(l_post)
}






compute_matrix<- function(n, sigma, K){
  Mat  = matrix(0,n,n )
  for (i in 1:n){
    for (j in 1:n){
      if (j<=i){
        julia_assign("i", as.integer(i))
        julia_assign("j", as.integer(j))
        julia_assign("H", as.integer(K))
        julia_assign("σ", sigma)
        Mat[i,j]= julia_eval( " log(GibbsTypePriors.Cnk(i, j, σ)) - j* log(H)|> Float64")
      }
    }
  }
  return (Mat)
}

Cnk<- function(j,v,H,sigma, n_k){
  julia_assign("n_k", as.integer(n_k))
  julia_assign("j", as.integer(j))
  julia_assign("σ", sigma)
  Cnk = julia_eval( " GibbsTypePriors.Cnk(n_k, j, σ) |> Float64")
   return(((v/H)^j)*Cnk)
}

#julia_assign("σ", sigma)
#Mat[i,j]= julia_eval( " log(GibbsTypePriors.Cnk(i, j, σ)) - j* log(H)|> Float64")

#c(15000, 1000, 200, 3800)

f_v <- function(v,H =10, n_vec = c(10, 15, 15, 6, 4),alpha=1,  sigma =0.5 ){
  f_vec<- c()
  K= length(n_vec)
  for (i in 1:K){
  #res= 0
  #  for (j in 1:n_vec[i]){
  #      res= res + ((v/H)^j)*Cnk(sigma, n_vec[i],j)
  #  }
    res_vec<- lapply(1:n_vec[i],Cnk,v=v, H=H, sigma=sigma, n_k=n_vec[i])
    f_vec[i ] = sum(unlist(res_vec))
  }
 return(exp(-v)* v^(alpha/sigma - 1) *prod(f_vec))
}


log_f_v <- function(v,H, n_vec,alpha,  sigma ){
  f_vec<- c()
  K= length(n_vec)
  for (i in 1:K){
    res= 0
    for (j in 1:n_vec[i]){
      julia_assign("n_k", as.integer(n_vec[i]))
      julia_assign("j", as.integer(j))
      julia_assign("σ", sigma)
      Cnk = julia_eval( " GibbsTypePriors.Cnk(n_k, j, σ) |> Float64")
      res= res + ((v/H)^j)*Cnk
      print(((v/H)^j)*Cnk)
    }
    f_vec[i ] = res
  }
  return(-v + (alpha/sigma - 1)*log(v) +  sum(log(f_vec)))
}


#n_v =c(15000, 1000, 200, 3800)

n_v =c(10, 15, 15, 6, 4)
alpha = 1
sigma = 0.5
H= 10

v_seq= seq(0,40, length.out = 30 )
#n_v =c(150, 10, 20, 38)


f_val <- lapply(v_seq, f_v, H=10, n_vec = n_v,alpha =1,  sigma = 0.5)
plot(v_seq, unlist(f_val))
plot(v_seq, ((v_seq)^2)*(unlist(f_val)))

#plot(v_seq, log(unlist(f_val)))

ptr_logv_mat <- create_xptr("log_v_pdf_comp_mat")
ptr_logv <- create_xptr("log_v_pdf_C")


Cnk_mat = compute_matrix(50, sigma = 0.5, K=10)
v_l = c()
for (i in 1:3000){
  v= ru_rcpp(logf = ptr_logv_mat,alpha=1, sigma=sigma,H=H,k = length(n_v), nk_vec=n_v,Cnk_mat=Cnk_mat, n=1,  d=1, init=1)
  v_l[i] = unlist(v$sim_vals[1])
}
for (i in 1:1000){
  v= ru_rcpp(logf = ptr_logv,alpha=1, sigma=sigma,H=H,k = length(n_v), nk_vec=n_v, n=1,  d=1, init=1)
  v_l[i] = unlist(v$sim_vals[1])
}

a = integrate(Vectorize(f_v), 0, Inf)

plot(v_seq, unlist(f_val)/a$value,col="green")
lines(density(v_l),type="l",col="red")

plot(density(v_l),type="l",col="red")







#optimize(f_v, lower = 0, upper = 10)


##### manual ratio of uniforms for logV








