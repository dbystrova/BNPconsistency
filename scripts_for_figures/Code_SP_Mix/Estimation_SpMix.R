############################################################################
####Gibbs Sampling for Fitting a Sparse Finite Multivariate Gaussian Mixture
############################################################################
library(inline)
library(rbenchmark)
library(frailtySurv)
library(rust)
library(RcppArmadillo)
library(arm)
library(Rcpp)
library(gtools)
library(JuliaCall)

Rcpp::sourceCpp('~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/src/user_fns.cpp')
#based on the blogpost "https://www.kent.ac.uk/smsas/personal/msr/rlaptrans.html"
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/rlaptrans.r")

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





sample.pym_w <- function(S_k, alpha, sigma, H, mat, fun =ptr_logv_mat  ){
  n_k = table(S_k)
  print(n_k)
  #sample v by ptr_logv_mat
  v_s = ru_rcpp(logf = fun,alpha=alpha, sigma=sigma,H=H,k = length(n_k), nk_vec=n_k,Cnk_mat=mat, n=1,  d=1, init=1)
  #sample lk
  print(v_s$sim_vals[1])
  lk <- sample_lk_mat(n_k,v_s$sim_vals[1],sigma,H,mat)
  #lh[c(as.numeric(c(names(n_k))))]= lk
  W_h <- rep(0,H)
  P_h<-  rep(0,H)
  p_vec<- n_k - lk*sigma
  alpha_post<- alpha + sum(lk)*sigma 
  
  W_h<- rdirichlet(1,c(p_vec, alpha_post))
  ### R
  Uv<- rgamma(1,alpha_post/sigma,alpha_post/sigma)
  U<- (Uv)^(1/sigma)
  x.rlap <- rlaptrans(H, lt.temp_st_pdf, c=alpha_post/(sigma*H), sigma, k=U)
  R_h<- x.rlap /sum(x.rlap)
  P_h[c(as.numeric(c(names(n_k))))]<- W_h[1:length(n_k)] + W_h[length(n_k)+1]* R_h[1:length(n_k)]
  P_h[-c(as.numeric(c(names(n_k))))] <-  W_h[length(n_k)+1]* R_h[(length(n_k)+1):H]
  return(P_h)
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
#estGibbs_2 <- MultVar_NormMixt_Gibbs_IndPriorNormalgamma(y, S_0, mu_0, sigma_0, eta_0, e0, c0, C0_0, 
 #                                                        g0, G0, b0, B0, nu, B_0, M, burnin, c_proposal, priorOnE0 = priorOnE0, lambda = FALSE,seed =seed_+1, sigma_py =  sigma_py)



MultVar_NormMixt_Gibbs_IndPriorNormalgamma <- function(y, S_0, mu_0, sigma_0, eta_0, e0, c0, C0_0, 
                                                       g0, G0, b0, B0k, nu, lam_0, M, burnin, c_proposal, priorOnE0, lambda,seed = 1, sigma_py =0) {
  
  set.seed(seed)
  print(seed)
  K <- ncol(mu_0)  #number of components 
  N <- nrow(y)  #number of observations
  r <- ncol(y)  #number of dimensions
  ##change!!!!
  R <- apply(y, 2, function(x) diff(range(x)))
  
  ## initializing current values:
  eta_0 <- rep(1/K, K)
  sigma_0 <- array(0, dim = c(r, r, K))
  mu_j <- matrix(0, r, K)
  S_j <- rep(0, N)
  C0_j <- matrix(0, r, r)
  
  eta_j <- eta_0
  sigma_j <- sigma_0
  invsigma_j <- sigma_j
  det_invsigma_j <- rep(0, K)
  mu_j <- mu_0
  S_j <- S_0
  #B_j <- B_0
  B_j <- lam_0
  
  C0_j <- C0_0
  #invB0_j <- solve(B0)
  invB0_j <- solve(B0k)
  b0_j <- b0
  Nk_j <- tabulate(S_j, K)
  print(Nk_j)
  e0_p <- 0
  #ptr_N01 <- create_xptr("log_v_pdf_C")
  ptr_logv_mat <- create_xptr("log_v_pdf_comp_mat")
  
  if (sigma_py > 0 ){
    
    julia <- julia_setup()
    julia_library("GibbsTypePriors")
    julia_library("DataFrames")
    julia_library("DataFramesMeta")
    
    Cnk_mat = compute_matrix(N, sigma_py, K )
 
  }
  print(Cnk_mat)
  ## generating matrices for storing the draws:
  result <- list(Eta = matrix(0, M, K), Mu = array(0, dim = c(M, r, K)), Sigma = array(0, dim = c(M, 
                                                                                                  r, r, K)), S_alt_matrix = matrix(0L, M, N), Nk_matrix_alt = matrix(0L, M, K), Nk_view = matrix(0L, 
                                                                                                                                                                                                 (M + burnin)/100 + 1, K), B = matrix(0, M, r), e0_vector = rep(0, M), mixlik = rep(0, M), 
                 mixprior = rep(0, M), nonnormpost = rep(0, M), cdpost = rep(0, M), nonnormpost_mode_list = vector("list", 
                                                                                                                   K), mixlik_mode_list = vector("list", K))
  
  
  acc <- rep(FALSE, M)  #for storing the acceptance results in the MH step
  
  
  ## Initialising the storing matrices:
  result$Mu[1, , ] <- mu_0
  result$Eta[1, ] <- eta_0
  result$S_alt_matrix[1, ] <- S_0
  #result$B[1, ] <- B_0
  result$B[1, ] <- lam_0
  result$Nk_matrix_alt[1, ] <- Nk_j
  for (k in 1:K) {
    result$nonnormpost_mode_list[[k]] <- list(nonnormpost = -(10)^18)
    result$mixlik_mode_list[[k]] <- list(mixlik = -(10)^18)
  }
  
  
  
  ### constant parameters for every iteration:
  p_gig <- nu - K/2
  a_gig <- 2 * nu
  gn <- g0 + K * c0
  
  
  
  
  ####################################################################################### simulation starts:
  s <- 1
  result$Nk_view[1, ] <- Nk_j
  m <- 2
  while (m <= M | m <= burnin) {
    
    if (m == burnin) {
      m <- 1
      burnin <- 0
    }
    
    
    ### relevant component specific quantities:
    mean_yk <- matrix(0, r, K)
    mean_yk <- sapply(1:K, function(k) colMeans(y[S_j == k, , drop = FALSE]))
    if (sum(is.na(mean_yk)) > 0) {
      ## to catch the case if a group is empty: NA values are substituted by zeros
      mean_yk[is.na(mean_yk)] <- 0
    }
    Nk_j <- tabulate(S_j, K)
    if (!(m%%100)) {
      cat("\n", m, " ", Nk_j)
      s <- s + 1
      result$Nk_view[s, ] <- Nk_j
    }
    Nk_alt_j <- Nk_j
    S_alt_j <- S_j
    K0_j <- sum(Nk_j != 0)  ##number of nonempty components
    
    
    #################### first step: parameter simulation (conditional on classification S_j): (1a): Sample eta_j:
    ek <- e0 + Nk_j
    if (sigma_py == 0){
      eta_j <- bayesm::rdirichlet(ek)
    } else{
      eta_j <- sample.pym_w(S_j, alpha = e0, sigma = sigma_py, H = K, mat =Cnk_mat, fun =ptr_logv_mat )
    }
    
    #### (1b): sample Sigma^{-1} for each component k: calculate posterior moments ck and Ck and sample
    #### from the inverted Wishart distribution:
    Ck <- array(0, dim = c(r, r, K))
    ck <- c0 + Nk_j/2
    for (k in 1:K) {
      if (Nk_j[k] != 0) {
        # Bettina:
        Ck[, , k] <- C0_j + 0.5 * crossprod(sweep(y[S_j == k, , drop = FALSE], 2, mu_j[, 
                                                                                       k], FUN = "-"))
      } else {
        Ck[, , k] <- C0_j
      }
      sig <- rwishart(2 * ck[k], 0.5 * chol2inv(chol(Ck[, , k])))  #attention: rwishart(nu,v)(Rossi)=> nu=2*c0,v=0.5*C0, wishart(c0,C0) (FS)
      sigma_j[, , k] <- sig$IW
      invsigma_j[, , k] <- sig$W
      det_invsigma_j[k] <- det(invsigma_j[, , k])
    }
    
    
    #### (1c): Sample mu_j for each component k:
    Bk <- array(0, dim = c(r, r, K))
    bk <- matrix(0, r, K)
    invB0_j <- diag(1/((R^2) * B_j))
    for (k in 1:K) {
      Bk[, , k] <- chol2inv(chol(invB0_j + invsigma_j[, , k] * Nk_j[k]))
      bk[, k] <- Bk[, , k] %*% (invB0_j %*% b0_j + invsigma_j[, , k] %*% mean_yk[, k] * Nk_j[k])
      mu_j[, k] <- t(chol(Bk[, , k])) %*% rnorm(r) + bk[, k]
    }
    
    
    
    #################### second step: classification of observations (conditional on knowing the parameters) Bettina:
    mat <- sapply(1:K, function(k) eta_j[k] * dmvnorm(y, mu_j[, k], sigma_j[, , k]))
    S_j <- apply(mat, 1, function(x) sample(1:K, 1, prob = x))
    Nk_j <- tabulate(S_j, K)
    Nk_neu_j <- Nk_j
    
    
    
    ###################### third step: sample the hyperparameters (3a): sample the hyperparameter-vector B conditionally
    ###################### on mu_j[,k] and sigma_j[,,k]:
    b_gig <- vector(length = r)
    b_gig <- (rowSums((mu_j - b0_j)^2))/R^2
    if (lambda == TRUE) {
      for (l in 1:r) {
        # accept/reject algorithm from Helga:
        if (b_gig[l] < 1e-06) {
          check <- 0
          while (check == 0) {
            ran <- 1/rgamma(1, shape = -p_gig, scale = 2/b_gig[l])
            check <- (runif(1) < exp(-a_gig/2 * ran))
          }
        } else {
          ### random generator from package 'Runuran': Create distribution object for GIG distribution
          distr <- udgig(theta = p_gig, psi = a_gig, chi = b_gig[l])  #b_gig=chi,a_gig=psi,p_gig=theta
          ## Generate generator object; use method PINV (inversion)
          gen <- pinvd.new(distr)
          ## Draw a sample of size 1
          ran <- ur(gen, 1)
        }
        B_j[l] <- ran
      }
    } else {
      #B_j <- B_0
      B_j <- lam_0
    }
    
    
    
    
    #### (3b): sample the hyperparameter C0 conditionally on sigma:
    C0_j <- rwishart(2 * gn, 0.5 * chol2inv(chol(G0 + rowSums(invsigma_j, dims = 2))))$W  #from package 'bayesm'
    
    
    
    #### (3c):assuming that the mean appearing in the normal prior on the group mean mu_k follows a
    #### improper prior p(b0)=const, sample b0 from N(1/K*sum(mu_i);1/K*B0 )
    B0 <- diag((R^2) * B_j)
    b0_j <- mvrnorm(1, rowSums(mu_j)/K, 1/K * B0)
    
    
    
    
    #### (3d): sample the hyperparameter e0 from p(e0|eta,a,b) via MH-step:
    a_gam <- 10
    b_gam <- a_gam  #e0~G(a_gam,b_gam*Kmax)
    # b_gam=1/K
    Kmax <- K
    const <- Kmax * b_gam
    ## current value:
    e0_j <- e0
    if (priorOnE0 == T) {
      ## proposal value:
      le0_p <- log(e0_j) + rnorm(1, 0, c_proposal)
      e0_p <- exp(le0_p)
      # likelihood ratio between the proposed value and the previous value:
      eta_j[eta_j == 0] <- 10^(-50)
      lalpha1 <- (e0_p - e0_j) * sum(log(eta_j)) + lgamma(K * e0_p) - lgamma(K * e0_j) - K * 
        (lgamma(e0_p) - lgamma(e0_j)) + (a_gam - 1) * (log(e0_p) - log(e0_j)) - (e0_p - e0_j) * 
        const + (log(e0_p) - log(e0_j))
      alpha1 <- min(exp(lalpha1), 1)
      ## 
      alpha2 <- runif(1)
      ## the proposed value is accepted with probability alpha2
      if (alpha2 <= alpha1) {
        e0_j <- e0_p
        if (burnin == 0) {
          acc[m] <- TRUE
        }
      }
    }
    result$e0_vector[m] <- e0_j
    e0 <- e0_j
    
    
    #### additional step (3e): evaluating the mixture likelihood and the complete-data posterior
    
    ## evaluating the mixture likelihood:
    mat_neu <- sapply(1:K, function(k) eta_j[k] * dmvnorm(y, mu_j[, k], sigma_j[, , k]))
    mixlik_j <- sum(log(rowSums(mat_neu)))
    
    ## evaluating the mixture prior:
    mixprior_j <- log(MCMCpack::ddirichlet(as.vector(eta_j), rep(e0, K))) + sum(dmvnorm(t(mu_j), 
                                                                                        b0, diag((R^2) * B_j), log = TRUE)) + sum(sapply(1:K, function(k) lndIWishart(2 * c0, 
                                                                                                                                                                      0.5 * C0_j, sigma_j[, , k]))) + lndIWishart(2 * g0, 0.5 * G0, C0_j) + dgamma(e0, shape = a_gam, 
                                                                                                                                                                                                                                                   scale = 1/b_gam, log = TRUE) + sum(dgamma(B_j, shape = nu, scale = 1/nu, log = TRUE))
    
    ## evaluating the nonnormalized complete-data posterior:
    cd <- c()
    for (k in 1:K) {
      if (sum(S_j == k)) {
        cd[S_j == k] <- dmvnorm(y[S_j == k, ], mu_j[, k], sigma_j[, , k], log = TRUE)
      }
    }
    cdpost_j <- sum(cd) + mixprior_j
    if (burnin == 0) {
      result$mixlik[m] <- mixlik_j
      result$mixprior[m] <- mixprior_j
      result$nonnormpost[m] <- result$mixlik[m] + result$mixprior[m]
      result$cdpost[m] <- cdpost_j
    }
    
    
    
    ###################### fourth step: permutation of the labeling and storing the results
    perm <- sample(K)
    
    if (burnin == 0) {
      ### storing the new values:
      result$Mu[m, , perm] <- mu_j
      result$Eta[m, perm] <- eta_j
      result$S_alt_matrix[m, ] <- perm[S_alt_j]
      result$Nk_matrix_alt[m, perm] <- Nk_alt_j
      result$B[m, ] <- B_j
      result$Sigma[m, , , perm] <- sigma_j
      result$acc_rate <- sum(acc)/M
    }
    
    if ((burnin == 0) & (result$nonnormpost[m] > result$nonnormpost_mode_list[[K0_j]]$nonnormpost)) {
      result$nonnormpost_mode_list[[K0_j]] <- list(nonnormpost = result$nonnormpost[m], mu = mu_j[, 
                                                                                                  Nk_alt_j != 0], bk = bk[, Nk_alt_j != 0], Bk = Bk[, , Nk_alt_j != 0], eta = eta_j[Nk_alt_j != 
                                                                                                                                                                                      0])
    }
    if ((burnin == 0) & (result$mixlik[m] > result$mixlik_mode_list[[K0_j]]$mixlik)) {
      result$mixlik_mode_list[[K0_j]] <- list(mixlik = result$mixlik[m], mu = mu_j[, Nk_alt_j != 
                                                                                     0], Sigma = sigma_j[, , Nk_alt_j != 0], eta = eta_j[Nk_alt_j != 0])  ##Bettina fragen:mixture likelihood oder complete data likelihood nehmen?
    }
    
    m <- m + 1
  }
  
  return(result)
}



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

