rm(list=ls()) 
graphics.off()
setwd("~/Documents/GitHub/BNPconsistency/data")

#### Try on BNPmix

library(BNPmix)
library(tidyverse)
library(latex2exp)
library(viridis)
library(MCMCpack)
library(JuliaCall)
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")

library(BSDA)
library(mclust)
data("thyroid")
data("Slc")
# data("diabetes")
# data("Galaxie")

dmvnorm <- function(x, mean, sigma, log = FALSE) {
  if (is.vector(x)) {
    x <- matrix(x, ncol = length(x))
  }
  if (missing(mean)) {
    mean <- rep(0, length = ncol(x))
  }
  if (missing(sigma)) {
    sigma <- diag(ncol(x))
  }
  if (NCOL(x) != NCOL(sigma)) {
    print(NCOL(x))
    print(NCOL(sigma))
    stop("x and sigma have non-conforming size")
  }
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (length(mean) != NROW(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  distval <- mahalanobis(x, center = mean, cov = chol2inv(chol(sigma)), inverted = TRUE)
  logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
  logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
  if (log) 
    return(logretval)
  exp(logretval)
}

compute_log_lik<- function(y, M, nburn, Pk, Mu, Sigma, univ){
  if(univ)
    n_data=length(y)
  else
    n_data = dim(y)[1]
  log_lik = c()
  for (i in 1:(M-nburn)){
    pk = Pk[[i]]; mu = Mu[[i]]; sig = Sigma[[i]]
    K = length(pk)
    lik = matrix(0,n_data,K )
    for (j in 1:K){
      if(univ)
        lik[,j] = pk[j]*dnorm(y,mu[j], sqrt(sig[j]))
      else
        lik[,j] = pk[j]*dmvnorm(y,mu[j,], sig[,,j])
    }
    log_lik[i] = log_lik[i] = sum(log(apply(lik,1,sum)))
  }
  return(log_lik)
}

l2norm_univ = function(theta.j, theta.i){
  l2_mu = sum((theta.j[[1]] -theta.i[[1]])^2)
  l2_Sigma = sum((theta.j[[2]] -theta.i[[2]])^2)
  return(sqrt(l2_mu + l2_Sigma))
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

apply_MTM_BNPmix<- function(fit_obj,c,n,  w_n = (log(n)/n)^(1/4), it = 10000){
  post_k<-list()
  it <- dim(fit_obj$clust)[1]
  n= dim(fit_obj$clust)[2]
  for (k in seq(2, it, 4)){
    G <- list()
    G$p = fit_obj$probs[[k]]
    G$theta = list()
    for (i in 1:length(fit_obj$probs[[k]])){
      #G$theta[[i]]<-  list(fit_obj$mean[[k]][i,],fit_obj$sigma2[[k]][i,]) #case for diagonal
      G$theta[[i]]<-  list(fit_obj$mean[[k]][i,],fit_obj$sigma2[[k]][,,i]) 
    }
    if(length(G$p) == length(G$theta))
      G.post <- MTM(G, w_n, c)
    else
      cat("iteration :", k, ", weights : ", length(G$p), ", atoms : ", length(G$theta), "\n")
    post_k <- append(post_k, length(G.post$p))
  }
  
  return(post_k)
}

apply_MTM_BNPmix_univ <- function(fit_obj, c, n,  w_n = (log(n)/n)^(1/4)){
  post_k<-list()
  it <- dim(fit_obj$clust)[1]
  n= dim(fit_obj$clust)[2]
  for (k in seq(2, it, 4)){
    G <- list()
    G$p = c(fit_obj$probs[[k]])
    # print(sum(G$p))
    G$theta = list()
    for (i in 1:length(fit_obj$mean[[k]])){
      G$theta[[i]]<-  list(fit_obj$mean[[k]][i,],fit_obj$sigma2[[k]][i,]) 
    }
    if(length(G$p) == length(G$theta))
      G.post <- MTM_univ(G, w_n, c)
    else
      cat("iteration :", k, ", weights : ", length(G$p), ", atoms : ", length(G$theta), "\n")
    post_k <- append(post_k, length(G.post$p))
  }
  
  return(post_k)
}

### Posterior for different alpha/sigma

py_mix <- function(it, burn, model = "LS", data, alpha, sigma, c_vec){
  mcmc <- list(niter = it, nburn = burn, model = model, method = "MAR")
  # Fit the PY mixture
  prior <- list(strength = alpha, discount = sigma, hyper = FALSE)
  output <- list(out_type = "FULL", out_param=TRUE)
  
  fit.sim1 <- PYdensity(y = data, mcmc = mcmc, prior = prior, output = output)
  fit.sim2 <- PYdensity(y = data, mcmc = mcmc, prior = prior, output = output)
  fit.sim3 <- PYdensity(y = data, mcmc = mcmc, prior = prior, output = output)
  fit.sim4 <- PYdensity(y = data, mcmc = mcmc, prior = prior, output = output)
  ## Visual diagnosis
  ll1<- compute_log_lik(data, it, burn, fit.sim1$probs, fit.sim1$mean, fit.sim1$sigma2, fit.sim1$univariate)
  ll2<- compute_log_lik(data, it, burn, fit.sim2$probs, fit.sim2$mean, fit.sim2$sigma2, fit.sim2$univariate)
  ll3<- compute_log_lik(data, it, burn, fit.sim3$probs, fit.sim3$mean, fit.sim3$sigma2, fit.sim3$univariate)
  ll4<- compute_log_lik(data, it, burn, fit.sim4$probs, fit.sim4$mean, fit.sim4$sigma2, fit.sim4$univariate)
  ## covergence diagnostics
  log_lik_combines <- mcmc.list(mcmc(ll1),mcmc(ll2),mcmc(ll3),mcmc(ll4))
  Rhat_ll<- gelman.diag(log_lik_combines)$psrf[1]
  print(paste("Rhat:", Rhat_ll))
  plot(c(ll1, ll2, ll3, ll4), type = "l", main = paste("Likelihood trace plot: alpha =", alpha, "and sig =", sigma))
  
  fun <-function(x){
    length(unique(x))
  } 
  clust_list1 = apply(fit.sim1$clust,1,fun)
  clust_list2 = apply(fit.sim2$clust,1,fun)
  clust_list = append(clust_list1,clust_list2)
  df= tibble(K= 1:max(clust_list), pk = rep(0,max(clust_list)))
  s= tibble(p_k = clust_list) %>%count(p_k) 
  df$pk[s$p_k] = s$n
  
  mean = rep(0,length(c_vec))
  map = rep(0,length(c_vec))
  if(is.null(dim(data))){
    n = length(data)
    for(i in 1:length(c_vec)){
      result_MTM1 = unlist(apply_MTM_BNPmix_univ(fit.sim1, c = c_vec[i], w_n = (log(n))^(-2)))
      result_MTM2 = unlist(apply_MTM_BNPmix_univ(fit.sim2, c = c_vec[i], w_n = (log(n))^(-2)))
      result_MTM3 = unlist(apply_MTM_BNPmix_univ(fit.sim3, c = c_vec[i], w_n = (log(n))^(-2)))
      result_MTM4 = unlist(apply_MTM_BNPmix_univ(fit.sim4, c = c_vec[i], w_n = (log(n))^(-2)))
      result_MTM = c(result_MTM1, result_MTM2, result_MTM3, result_MTM4)
      map[i] = as.numeric(names(table(result_MTM)[which.max(table(result_MTM))]))
      mean[i] = mean(round(result_MTM))
    }
  }
  else{
    n = dim(data)[1]
    for(i in 1:length(c_vec)){
      result_MTM1 = unlist(apply_MTM_BNPmix(fit.sim1, c = c_vec[i], w_n = (log(n))^(-2)))
      result_MTM2 = unlist(apply_MTM_BNPmix(fit.sim2, c = c_vec[i], w_n = (log(n))^(-2)))
      result_MTM3 = unlist(apply_MTM_BNPmix(fit.sim3, c = c_vec[i], w_n = (log(n))^(-2)))
      result_MTM4 = unlist(apply_MTM_BNPmix(fit.sim4, c = c_vec[i], w_n = (log(n))^(-2)))
      result_MTM = c(result_MTM1, result_MTM2, result_MTM3, result_MTM4)
      map[i] = as.numeric(names(table(result_MTM)[which.max(table(result_MTM))]))
      mean[i] = mean(round(result_MTM))
    }
  }
  MTM_out = tibble(c=c_vec, mean = mean, map = map)
  return(list(df=df, MTM_out=MTM_out))
}


###### Univariate PY #######
alpha_seq = c(0.01,0.5)
sigma_seq = c(0.1,0.25)

it = 20000; burn = 10000
c_list = list(seq(0, 3, length.out = 50), seq(0, 5, length.out = 50), seq(0, 5.5, length.out = 50), seq(0, 6, length.out = 50))
data= unname(as_vector(Slc))
n = length(data)

for (j in 1:length(alpha_seq)){
  for (k in 1:length(sigma_seq)){
    c_vec = c_list[[k-1+(2*j-1)]]
    if ((alpha_seq[j] == alpha_seq[1])&&(sigma_seq[k] == sigma_seq[1])){
      print(paste("alpha : ", alpha_seq[j], " sigma : ", sigma_seq[k]))
      out = py_mix(it =it, burn = burn, data= data, alpha = alpha_seq[j], sigma = sigma_seq[k], c_vec=c_vec)
      df = out$df
      df$alpha = rep(alpha_seq[j], dim(df)[1])
      df$sigma = rep(sigma_seq[k], dim(df)[1])
      df$Process_type = paste0("Pkn_",k-1+(2*j-1))
      MTM_out = out$MTM_out
      MTM_out$alpha = rep(alpha_seq[j], dim(MTM_out)[1])
      MTM_out$sigma = rep(sigma_seq[k], dim(MTM_out)[1])
      MTM_out$Process_type = paste0("Pkn_",k-1+(2*j-1))
    }
    else {
      print(paste("alpha : ", alpha_seq[j], " sigma : ", sigma_seq[k]))
      out = py_mix(it =it, burn = burn, data= data, alpha = alpha_seq[j], sigma = sigma_seq[k], c_vec=c_vec)
      df_ = out$df
      df_$alpha = rep(alpha_seq[j], dim(df_)[1])
      df_$sigma = rep(sigma_seq[k], dim(df_)[1])
      df_$Process_type = paste0("Pkn_",k-1+(2*j-1))
      df =rbind(df, df_)
      MTM_out_ = out$MTM_out
      MTM_out_$alpha = rep(alpha_seq[j], dim(MTM_out_)[1])
      MTM_out_$sigma = rep(sigma_seq[k], dim(MTM_out_)[1])
      MTM_out_$Process_type = paste0("Pkn_",k-1+(2*j-1))
      MTM_out =rbind(MTM_out, MTM_out_)
    }
  }
}

names(MTM_out) <- c('c', 'Mean', 'MAP', 'alpha', 'sigma', 'Process_type')
MTM_out <- MTM_out%>%gather(key = type, value = val, Mean, MAP)%>%
  group_by(Process_type,alpha,sigma)

df <- df%>%group_by(Process_type,alpha,sigma)%>%mutate(pk =pk/sum(pk))
julia <- julia_setup()
julia_library("GibbsTypePriors")

K_= max(df$K)
julia_assign("K_bound", K_)
julia_assign("A1", alpha_seq[1])
julia_assign("A2", alpha_seq[2])
julia_assign("S1", sigma_seq[1])
julia_assign("S2", sigma_seq[2])
# alpha = e0*K
julia_assign("N", n)
# print( fig_df$line$alpha[1]*K_)
#a = julia_eval("Pkn_Dirichlet_mult.(1:K_bound,N,K_bound, A1)")

df_prior = tibble(K= 1:K_, 
                  Pkn_1 = round(julia_eval("Pkn_2PD.(1:K_bound, N, A1, S1)"),3),
                  Pkn_2 = round(julia_eval("Pkn_2PD.(1:K_bound, N, A1, S2)"),3), 
                  Pkn_3 = round(julia_eval("Pkn_2PD.(1:K_bound, N, A2, S1)"),3),
                  Pkn_4 = round(julia_eval("Pkn_2PD.(1:K_bound, N, A2, S2)"),3))%>% gather(Process_type, pk,Pkn_1:Pkn_4)

df_prior$Type =rep("Prior", dim(df_prior)[1])
df$Type= rep("Posterior", dim(df)[1])

df_merged = rbind(df_prior,df[, c("K","Process_type", "pk", "Type")] )

ph <- ggplot(df_merged, aes(K,pk,color =Process_type))+geom_bar(aes(linetype=Type),size = 0.7, stat="identity",alpha =0.0, position = "identity")+
  geom_vline(xintercept=2,  linetype="dashed")+ylab("Density")+xlab(TeX('$K_n$'))+
  scale_y_continuous(limits = c(0,1))+
  scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D",name = TeX(sprintf('$\\alpha, \\sigma$')), labels=unname(TeX(c(sprintf('$\\alpha$ = %.2f, $\\sigma$ = %.1f', alpha_seq[1], sigma_seq[1]),
                                                                                                                                         sprintf('$\\alpha$ = %.2f, $\\sigma$ = %.2f', alpha_seq[1], sigma_seq[2]),
                                                                                                                                         sprintf('$\\alpha$ = %.1f, $\\sigma$ = %.1f', alpha_seq[2], sigma_seq[1]),
                                                                                                                                         sprintf('$\\alpha$ = %.1f, $\\sigma$ = %.2f', alpha_seq[2], sigma_seq[2])))))+
  scale_linetype_manual(name = "Distribution",values = c(1, 2),guide = guide_legend(override.aes = list(linetype = c(1, 2),color = "black") ) )+
  theme_minimal()+
  theme(strip.background = element_blank(), strip.text.x = element_blank())+ facet_wrap(~Process_type)
plot(ph)
# pdf(paste0("Figure_data_alpha_sigma",fig_df_mut$alpha[1],"_3_.pdf" ))
pdf(paste0("Figure_Slc_alpha_sigma_grid_L.pdf" ))
plot(ph)
dev.off()



pm <- ggplot(MTM_out, aes(x = c, y = val, color = Process_type))+ geom_line(aes(linetype=type))+
  scale_linetype_manual(name = "Distribution",values = c("solid", "dotdash"),guide = guide_legend(override.aes = list(color = "black") ) )+
  scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D",name = TeX(sprintf('$\\alpha, \\sigma$')), labels=unname(TeX(c(sprintf('$\\alpha$ = %.2f, $\\sigma$ = %.1f', alpha_seq[1], sigma_seq[1]),
                                                                                                                                         sprintf('$\\alpha$ = %.2f, $\\sigma$ = %.2f', alpha_seq[1], sigma_seq[2]),
                                                                                                                                         sprintf('$\\alpha$ = %.1f, $\\sigma$ = %.1f', alpha_seq[2], sigma_seq[1]),
                                                                                                                                         sprintf('$\\alpha$ = %.1f, $\\sigma$ = %.2f', alpha_seq[2], sigma_seq[2])))))+
  
  theme_minimal()+  theme(strip.background = element_blank(), strip.text.x = element_blank())+
  ylab(expression(tilde(K)))+xlab(TeX('$c$')) +geom_hline(yintercept = 2, linetype="dashed", size = 0.5, alpha =0.5) +    facet_wrap(~Process_type)
plot(pm)
pdf(paste0("Figure_Slc_PY_MTM.pdf"))
plot(pm)
dev.off()

###### Multivariate PY #######
it = 70000; burn = 50000
c_vec = seq(0, 13, length.out = 50)
data= thyroid[2:6]
n = dim(data)[1]

# fit.sim_list = list()
for (j in 1:length(alpha_seq)){
  for (k in 1:length(sigma_seq)){
    if ((alpha_seq[j] == alpha_seq[1])&&(sigma_seq[k] == sigma_seq[1])){
      print(paste("alpha : ", alpha_seq[j], " sigma : ", sigma_seq[k]))
      out = py_mix(it =it, burn = burn, data= data, alpha = alpha_seq[j], sigma = sigma_seq[k], c_vec=c_vec)
      df = out$df
      df$alpha = rep(alpha_seq[j], dim(df)[1])
      df$sigma = rep(sigma_seq[k], dim(df)[1])
      df$Process_type = paste0("Pkn_",k-1+(2*j-1))
      MTM_out = out$MTM_out
      MTM_out$alpha = rep(alpha_seq[j], dim(MTM_out)[1])
      MTM_out$sigma = rep(sigma_seq[k], dim(MTM_out)[1])
      MTM_out$Process_type = paste0("Pkn_",k-1+(2*j-1))
      # fit.sim_list[[j*length(sigma_seq)-1+k-1]]=out$fit.sim
    }
    else {
      print(paste("alpha : ", alpha_seq[j], " sigma : ", sigma_seq[k]))
      out = py_mix(it =it, burn = burn, data= data, alpha = alpha_seq[j], sigma = sigma_seq[k], c_vec=c_vec)
      df_ = out$df
      df_$alpha = rep(alpha_seq[j], dim(df_)[1])
      df_$sigma = rep(sigma_seq[k], dim(df_)[1])
      df_$Process_type = paste0("Pkn_",k-1+(2*j-1))
      df =rbind(df, df_)
      MTM_out_ = out$MTM_out
      MTM_out_$alpha = rep(alpha_seq[j], dim(MTM_out_)[1])
      MTM_out_$sigma = rep(sigma_seq[k], dim(MTM_out_)[1])
      MTM_out_$Process_type = paste0("Pkn_",k-1+(2*j-1))
      MTM_out =rbind(MTM_out, MTM_out_)
      # fit.sim_list[[j*length(sigma_seq)-1+k-1]]=out$fit.sim
    }
  }
}

names(MTM_out) <- c('c', 'Mean', 'MAP', 'alpha', 'sigma', 'Process_type')
MTM_out <- MTM_out%>%gather(key = type, value = val, Mean, MAP)%>%
  group_by(Process_type,alpha,sigma)

df <- df%>%group_by(Process_type,alpha,sigma)%>%mutate(pk =pk/sum(pk))
julia <- julia_setup()
julia_library("GibbsTypePriors")

K_= max(df$K)
julia_assign("K_bound", K_)
julia_assign("A1", alpha_seq[1])
julia_assign("A2", alpha_seq[2])
julia_assign("S1", sigma_seq[1])
julia_assign("S2", sigma_seq[2])
# alpha = e0*K
julia_assign("N", n)
# print( fig_df$line$alpha[1]*K_)
#a = julia_eval("Pkn_Dirichlet_mult.(1:K_bound,N,K_bound, A1)")

df_prior = tibble(K= 1:K_, 
                  Pkn_1 = round(julia_eval("Pkn_2PD.(1:K_bound, N, A1, S1)"),3),
                  Pkn_2 = round(julia_eval("Pkn_2PD.(1:K_bound, N, A1, S2)"),3), 
                  Pkn_3 = round(julia_eval("Pkn_2PD.(1:K_bound, N, A2, S1)"),3),
                  Pkn_4 = round(julia_eval("Pkn_2PD.(1:K_bound, N, A2, S2)"),3))%>% gather(Process_type, pk,Pkn_1:Pkn_4)

df_prior$Type =rep("Prior", dim(df_prior)[1])
df$Type= rep("Posterior", dim(df)[1])

df_merged = rbind(df_prior,df[, c("K","Process_type", "pk", "Type")] )

ph <- ggplot(df_merged, aes(K,pk,color =Process_type))+geom_bar(aes(linetype=Type),size = 0.7, stat="identity",alpha =0.0, position = "identity")+
  geom_vline(xintercept=3,  linetype="dashed")+ylab("Density")+xlab(TeX('$K_n$'))+
  scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D",name = TeX(sprintf('$\\alpha, \\sigma$')), labels=unname(TeX(c(sprintf('$\\alpha$ = %.2f, $\\sigma$ = %.1f', alpha_seq[1], sigma_seq[1]),
                                                                                                                                         sprintf('$\\alpha$ = %.2f, $\\sigma$ = %.2f', alpha_seq[1], sigma_seq[2]),
                                                                                                                                         sprintf('$\\alpha$ = %.1f, $\\sigma$ = %.1f', alpha_seq[2], sigma_seq[1]),
                                                                                                                                         sprintf('$\\alpha$ = %.1f, $\\sigma$ = %.2f', alpha_seq[2], sigma_seq[2])))))+
  scale_linetype_manual(name = "Distribution",values = c(1, 2),guide = guide_legend(override.aes = list(linetype = c(1, 2),color = "black") ) )+
  theme_minimal()+
  theme(strip.background = element_blank(), strip.text.x = element_blank())+ facet_wrap(~Process_type)
plot(ph)
# pdf(paste0("Figure_data_alpha_sigma",fig_df_mut$alpha[1],"_3_.pdf" ))
pdf(paste0("Figure_thyroid_alpha_sigma_grid_L.pdf" ))
plot(ph)
dev.off()



pm <- ggplot(MTM_out, aes(x = c, y = val, color = Process_type))+ geom_line(aes(linetype=type))+scale_linetype_manual(values=c("solid", "dotdash"))+
  scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D",name = TeX(sprintf('$\\alpha, \\sigma$')), labels=unname(TeX(c(sprintf('$\\alpha$ = %.2f, $\\sigma$ = %.1f', alpha_seq[1], sigma_seq[1]),
                                                                                                                                         sprintf('$\\alpha$ = %.2f, $\\sigma$ = %.2f', alpha_seq[1], sigma_seq[2]),
                                                                                                                                         sprintf('$\\alpha$ = %.1f, $\\sigma$ = %.1f', alpha_seq[2], sigma_seq[1]),
                                                                                                                                         sprintf('$\\alpha$ = %.1f, $\\sigma$ = %.2f', alpha_seq[2], sigma_seq[2])))))+
  
  theme_minimal()+theme(strip.background = element_blank(), strip.text.x = element_blank())+ 
  ylab(expression(tilde(K)))+xlab(TeX('$c$')) +geom_hline(yintercept = 3, linetype="dashed", size = 0.5, alpha =0.5) +    facet_wrap(~Process_type)
plot(pm)
pdf(paste0("Figure_thyroid_PY_MTM.pdf"))
plot(pm)
dev.off()
