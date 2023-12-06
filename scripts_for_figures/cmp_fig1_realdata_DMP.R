source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Random_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Estimation_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Identification_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Gibbs_sampling_function.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/plt_Figure1.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")
library(JuliaCall)


comparison_data_alpha<- function(data,K_, M_it, nburn, alpha_l){
  pk <- list()
  R_h <- c()
  W_non_sorted <- list()
  W <- list()
  Mu_mat <- list()
  S_mat <- list()
  
  N <- dim(data$y)[1]
  for (i in 1:length(alpha_l)){
    pk[[i]] <- MCMC_function(data, e0=alpha_l[i], K=K_, M=M_it, burnin=nburn)  #conversion from Dir to DPM 
    Eta_ <- matrix(NA, nrow=dim(pk[[i]]$Eta)[1], ncol=dim(pk[[i]]$Eta)[2])
    for (j in 1:dim(pk[[i]]$Eta)[1]){
      Eta_[j,] <- sort(pk[[i]]$Eta[j,], decreasing=TRUE)
    }
    W[[i]] <- Eta_
    W_non_sorted[[i]] <- pk[[i]]$Eta
    Mu_mat[[i]] <- pk[[i]]$Mu
    S_mat[[i]] <-pk[[i]]$Sigma
  }
  df_ <- tibble(K=1:K_)
  df3_ <- tibble(K=1:K_)
  df4_ <- tibble(K=1:K_)
  for (j in 1:length(alpha_l)){
    name_ <- paste("Pkn_", j, sep="")
    name3_ <- paste("N_", j, sep="")
    name4_ <- paste("Alpha_", j, sep="")
    name5_ <- paste("W_", j, sep="")
    df_[,name_] <- pk[[j]]$p_k
    df3_[,name3_] <- rep(N, length(pk[[j]]$p_k))
    df4_[,name4_] <- rep(alpha_l[j]/K_, length(pk[[j]]$p_k))
  }
  df <- df_ %>% gather(Process_type, density, paste("Pkn_", 1, sep=""):paste("Pkn_", length(alpha_l), sep=""))
  df3 <- df3_ %>% gather(N_, N_val, paste("N_", 1, sep=""):paste("N_", length(alpha_l), sep=""))
  df4 <- df4_ %>% gather(Alpha_, Alpha_val, paste("Alpha_", 1, sep=""):paste("Alpha_", length(alpha_l), sep=""))
  
  df$Al <- df4$Alpha_val
  df$N <- df3$N_val
  
  df_l_ <- tibble(K= 1:((M_it)*2))
  df3_l_ <- tibble(K= 1:((M_it)*2))
  df4_l_ <- tibble(K= 1:((M_it)*2))
  for (j in 1:length(alpha_l)){
    name_ <- paste("P_", j, sep="")
    name3_ <- paste("N_", j, sep="")
    name4_ <- paste("Alpha_", j, sep="")
    df_l_[,name_] <- pk[[j]]$p_k_all
    df3_l_[,name3_] <- rep(N, length(pk[[j]]$p_k_all))
    df4_l_[,name4_] <- rep(alpha_l[j], length(pk[[j]]$p_k_all))
  }
  df_l <- df_l_ %>% gather(Process_type, density, paste("P_", 1, sep=""):paste("P_", length(alpha_l), sep=""))
  df_l3 <- df3_l_ %>% gather(N_, N_val, paste("N_", 1, sep=""):paste("N_", length(alpha_l), sep=""))
  df_l4 <- df4_l_ %>% gather(Alpha_, Alpha_val, paste("Alpha_", 1, sep=""):paste("Alpha_", length(alpha_l), sep=""))
  df_l$N <- df_l3$N_val
  df_l$Al <- df_l4$Alpha_val
  return(list(line = df, hist = df_l, eta = W_non_sorted, mu = Mu_mat, sig = S_mat))
}

alpha = 1
K_ = 10

##################################
######## Multivariate DMP ########
##################################
M_it = 10000
nburn = 30000

data_raw <- loadRData("~/Documents/GitHub/BNPconsistency/scripts_for_figures/real_data/thyroid.RData")
data <- list()
data$y <- data_raw[,2:6]
N <- dim(data_raw[,2:6])[1]

alpha_l <- c(0.01,0.5, 1, 10)
final <- comparison_data_alpha(data = data, K_ = 10 , M_it = M_it, nburn = nburn, alpha_l = alpha_l)
save(final, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_thyroid_DMP.RData")


################################
######## Univariate DMP ########
################################
M_it = 20000
nburn = 30000

data_raw <- loadRData("~/Documents/GitHub/BNPconsistency/scripts_for_figures/real_data/Slc.RData")
data <- list()
data$y <- unname(data_raw)
N <- dim(data_raw)[1]

alpha_l <- c(0.01, 0.5,1, 2)
final <- comparison_data_alpha(data = data, K_ = 10 , M_it = M_it, nburn = nburn, alpha_l = alpha_l)
save(final, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_Slc_DMP.RData")

