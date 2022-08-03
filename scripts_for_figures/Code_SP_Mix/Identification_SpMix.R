###########################################################################################################
## Clustering of the Mu-draws in the point process representation, with k-means based on Mahalanobis distance.
## The obtained labeling  is used toreorder the Eta draws and the allocation-matrix S.
## Only iterations resulting in 'true' permutations  are used for parameter estimation.
###########################################################################################################



MultVar_clust_FS_FAST_2 <- function(Mu, Eta, S_matrix, map_mu, map_cov, maha) {
  
  K <- length(Mu[1, 1, ])
  M <- length(Mu[, 1, 1])
  r <- length(Mu[1, , 1])
  N <- length(S_matrix[1, ])
  
  
  ############################## step 1:Cluster the MCMC draws create a matrix of size (M*K)xr, 
  # formed from all MCMC-draws of Mu, for applying a k-means clustering algorithm
  KM <- (Mu[, , 1])
  for (k in 2:K) {
    KM <- rbind(KM, Mu[, , k])
  }
  names_Mu <- c()
  for (j in 1:r) {
    names_Mu <- c(names_Mu, paste("mu", j, sep = ""))
  }
  colnames(KM) <- c(names_Mu)
  
  #### cluster Mu-draws with k-means, based on Mahalanobis-distance or without
  if (maha == TRUE) {
    #### 1a) flexclust-function are rewritten
    centMaha <- function(x) c(colMeans(x), as.vector(var(x)))
    distMaha <- function(x, centers) {
      z <- matrix(0, nrow = nrow(x), ncol = nrow(centers))
      for (k in 1:nrow(centers)) {
        mu <- centers[k, seq_len(ncol(x))]
        Sigma <- matrix(centers[k, -seq_len(ncol(x))], nrow = ncol(x))
        z[, k] <- sqrt(mahalanobis(x, mu, Sigma))
      }
      z
    }
    allcentMaha <- function(x, cluster, k = max(cluster, na.rm = TRUE)) {
      centers <- matrix(NA, nrow = k, ncol = ncol(x) + ncol(x)^2)
      for (n in 1:k) {
        if (sum(cluster == n, na.rm = TRUE) > 0) {
          centers[n, ] <- centMaha(x[cluster == n, , drop = FALSE])
        }
      }
      centers
    }
    kccafamily <- kccaFamily(cent = centMaha, dist = distMaha)
    kccafamily@allcent <- allcentMaha
    #### 1b) Mu-draws, arranged in one 'column', are picked up
    x <- KM
    colnames(x) <- c()
    #### 1c) staring classification is obtained by calculating the maha-distance between draws and
    # map-estimator of mu_k and cov_k
    center_start <- matrix(NA, nrow = k, ncol = ncol(x) + ncol(x)^2)
    for (n in 1:K) {
      center_start[n, ] <- c(map_mu[, n], as.vector(map_cov[, , n]))
    }
    class <- max.col(-distMaha(x, center_start))
    #### 1d) now the 'real'' classification is performed
    kmaha <- kcca(x, class, family = kccafamily)
    
    #### 1e) classification matrix is constructed:
    Rho_m <- c()
    for (l in 0:(K - 1)) {
      Rho_m <- cbind(Rho_m, clusters(kmaha)[(l * M + 1):((l + 1) * M)])
    }
    
  } else {
    #### standard clustering with kcca: 1a'') not necessary 1b'') Mu-draws, arranged in one 'column' are
    # picked up
    x <- KM
    colnames(x) <- c()
    #### 1c'') staring classification is obtained by k-means clustering with known cluster
    #### means(=mu_map)
    cent <- t(map_mu)
    cl_y <- kmeans(x, centers = cent, nstart = 30)
    class <- cl_y$cluster
    #### 1d'') now the (standard) classification
    kstand <- kcca(x, cent)
    par(mfrow = c(1, 1))
    #### 1e') classification matrix is constructed:
    Rho_m <- c()
    for (l in 0:(K - 1)) {
      Rho_m <- cbind(Rho_m, clusters(kstand)[(l * M + 1):((l + 1) * M)])
    }
  }
  
  ######################### step 2: Identifying non-permutations searching for permutations:
  m_rho <- c()
  for (m in 1:M) {
    if (any(sort(Rho_m[m, ]) != 1:K)) 
      m_rho <- c(m_rho, m)
  }
  non_perm_rate <- length(m_rho)/M
  non_perm_rate  ###rate of non-permutations
  
  
  ######################### step3: Relabel draws of Mu, Eta and S: unique labeling is achieved by reordering the draws
  # trough Rho_m:
  Mu_reord <- array(0, dim = c(M, r, K))
  Eta_reord <- matrix(0, M, K)
  S_matrix_reord <- matrix(0, M, N)
  
  for (m in 2:M) {
    Mu_reord[m, , Rho_m[m, ]] <- Mu[m, , ]
    Eta_reord[m, Rho_m[m, ]] <- Eta[m, ]
    S_matrix_reord[m, ] <- Rho_m[m, ][S_matrix[m, ]]  ##
  }
  
  ######################### step4:drop draws which are not permutations:
  only_2_permutations <- 0
  if (non_perm_rate < 0.999) {
    Mu_only_perm <- Mu_reord[setdiff(1:M, m_rho), , ]
    Eta_only_perm <- Eta_reord[setdiff(1:M, m_rho), ]
    S_matrix_only_perm <- S_matrix_reord[setdiff(1:M, m_rho), ]
  } else {
    Mu_only_perm <- Mu_reord
    Eta_only_perm <- Eta_reord
    S_matrix_only_perm <- S_matrix_reord
    only_2_permutations <- only_2_permutations + 1
  }
  
  return(list(S_matrix_only_perm = S_matrix_only_perm, Mu_only_perm = Mu_only_perm, Eta_only_perm = Eta_only_perm, 
              non_perm_rate = non_perm_rate, Rho_m = Rho_m, m_rho = m_rho))
}

