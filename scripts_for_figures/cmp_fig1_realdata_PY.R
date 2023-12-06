library(BNPmix)
library(tidyverse)
library(latex2exp)
library(viridis)
library(MCMCpack)
library(JuliaCall)
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")

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


### Posterior for different alpha/sigma

py_mix <- function(it, burn, model = "LS", data, alpha, sigma){
  mcmc <- list(niter = it, nburn = burn, model = model, method = "MAR")
  # Fit the PY mixture
  prior <- list(strength = alpha, discount = sigma, hyper = FALSE)
  output <- list(out_type = "FULL", out_param=TRUE)
  
  fit.sim1 <- PYdensity(y = data, mcmc = mcmc, prior = prior, output = output)
  fit.sim2 <- PYdensity(y = data, mcmc = mcmc, prior = prior, output = output)
  fit.sim3 <- PYdensity(y = data, mcmc = mcmc, prior = prior, output = output)
  fit.sim4 <- PYdensity(y = data, mcmc = mcmc, prior = prior, output = output)
  ## con vergence diagnosis
  ll1 <- compute_log_lik(data, it, burn, fit.sim1$probs, fit.sim1$mean, fit.sim1$sigma2, fit.sim1$univariate)
  ll2 <- compute_log_lik(data, it, burn, fit.sim2$probs, fit.sim2$mean, fit.sim2$sigma2, fit.sim2$univariate)
  ll3 <- compute_log_lik(data, it, burn, fit.sim3$probs, fit.sim3$mean, fit.sim3$sigma2, fit.sim3$univariate)
  ll4 <- compute_log_lik(data, it, burn, fit.sim4$probs, fit.sim4$mean, fit.sim4$sigma2, fit.sim4$univariate)
  log_lik_combines <- mcmc.list(mcmc(ll1),mcmc(ll2),mcmc(ll3),mcmc(ll4))
  Rhat_ll<- gelman.diag(log_lik_combines)$psrf[1]
  print(paste("Rhat:", Rhat_ll))
  plot(c(ll1, ll2, ll3, ll4), type = "l", main = paste("Likelihood trace plot: alpha =", alpha, "and sig =", sigma))
  
  fit.sim1 <- list(clust=fit.sim1$clust, probs=fit.sim1$probs, mean=fit.sim1$mean, sigma2=fit.sim1$sigma2)
  fit.sim2 <- list(clust=fit.sim2$clust, probs=fit.sim2$probs, mean=fit.sim2$mean, sigma2=fit.sim2$sigma2)
  fit.sim3 <- list(clust=fit.sim3$clust, probs=fit.sim3$probs, mean=fit.sim3$mean, sigma2=fit.sim3$sigma2)
  fit.sim4 <- list(clust=fit.sim4$clust, probs=fit.sim4$probs, mean=fit.sim4$mean, sigma2=fit.sim4$sigma2)
  
  if(is.null(dim(data))){
    n <- length(data)
    univ <- TRUE
  }
  else{
    n <- dim(data)[1]
    univ <- FALSE
  }
  return(list(f1 = fit.sim1, f2 = fit.sim2, f3 = fit.sim3, f4 = fit.sim4, Al = alpha, Sig = sigma, N = n, univ = univ))
}

alpha_seq = c(0.01,0.5)
sigma_seq = c(0.1,0.25)

#################################
######## Multivariate PY ########
#################################
it = 70000; burn = 60000

data_raw <- loadRData("~/Documents/GitHub/BNPconsistency/scripts_for_figures/real_data/thyroid.RData")
data <- data_raw[,2:6]

final <- list()
for (j in 1:length(alpha_seq)){
  for (k in 1:length(sigma_seq)){
      print(paste("alpha : ", alpha_seq[j], " sigma : ", sigma_seq[k]))
      out <- py_mix(it = it, burn = burn, data = data, alpha = alpha_seq[j], sigma = sigma_seq[k])
      final[[j*length(sigma_seq)-1+k-1]] <- out
  }
}

save(final, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_thyroid_PY.RData")

###############################
######## Univariate PY ########
###############################
it = 20000; burn = 10000

data_raw <- loadRData("~/Documents/GitHub/BNPconsistency/scripts_for_figures/real_data/Slc.RData")
data <- unname(as_vector(data_raw))

final <- list()
# fit.sim_list = list()
for (j in 1:length(alpha_seq)){
  for (k in 1:length(sigma_seq)){
    print(paste("alpha : ", alpha_seq[j], " sigma : ", sigma_seq[k]))
    out <- py_mix(it = it, burn = burn, data = data, alpha = alpha_seq[j], sigma = sigma_seq[k])
    final[[j*length(sigma_seq)-1+k-1]] <- out
  }
}

save(final, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_Slc_PY.RData")
