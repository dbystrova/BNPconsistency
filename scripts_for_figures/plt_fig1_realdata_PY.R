## read sources
library(tidyverse)
library(latex2exp)
library(viridis)
library(JuliaCall)
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")

fig_path <- "~/Documents/GitHub/BNPconsistency/figures/Figure_real_data/"

to_df <- function(final){
  for(i in 1:length(final)){
    fin_i <- final[[i]]
    clust_list <- apply(fin_i$f1$clust,1,function(x){length(unique(x))}) %>%
      append(apply(fin_i$f2$clust,1,function(x){length(unique(x))})) %>%
      append(apply(fin_i$f3$clust,1,function(x){length(unique(x))})) %>%
      append(apply(fin_i$f4$clust,1,function(x){length(unique(x))}))
    
    if(i == 1){
      df <- tibble(K = 1:max(clust_list), pk = rep(0, max(clust_list)))
      s <- tibble(p_k = clust_list) %>% count(p_k) 
      df$pk[s$p_k] <- s$n
      df$alpha <- rep(fin_i$Al, dim(df)[1])
      df$sigma <- rep(fin_i$Sig, dim(df)[1])
      df$Process_type <- paste0("Pkn_",i)
    }
    else{
      df_ <- tibble(K = 1:max(clust_list), pk = rep(0, max(clust_list)))
      s <- tibble(p_k = clust_list) %>% count(p_k) 
      df_$pk[s$p_k] <- s$n
      df_$alpha <- rep(fin_i$Al, dim(df_)[1])
      df_$sigma <- rep(fin_i$Sig, dim(df_)[1])
      df_$Process_type <- paste0("Pkn_",i)
      df <- rbind(df, df_)
    }
  }
  
  return(df)
}

######################################
######## Plot Multivariate PY ########
######################################
final <- loadRData( "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_thyroid_PY.RData")
df <- to_df(final)

alpha_seq <- as.numeric(levels(as.factor(df$alpha)))
sigma_seq <- as.numeric(levels(as.factor(df$sigma)))
N <- final[[1]]$N

df <- df%>%group_by(Process_type,alpha,sigma)%>%mutate(pk = pk/sum(pk))
julia <- julia_setup()
julia_library("GibbsTypePriors")

K_ <- max(df$K)
julia_assign("K_bound", K_)
julia_assign("A1", alpha_seq[1])
julia_assign("A2", alpha_seq[2])
julia_assign("S1", sigma_seq[1])
julia_assign("S2", sigma_seq[2])
julia_assign("N", N)

df_prior <- tibble(K = 1:K_, 
                  Pkn_1 = round(julia_eval("Pkn_2PD.(1:K_bound, N, A1, S1)"),3),
                  Pkn_2 = round(julia_eval("Pkn_2PD.(1:K_bound, N, A1, S2)"),3), 
                  Pkn_3 = round(julia_eval("Pkn_2PD.(1:K_bound, N, A2, S1)"),3),
                  Pkn_4 = round(julia_eval("Pkn_2PD.(1:K_bound, N, A2, S2)"),3)) %>% 
  gather(Process_type, pk,Pkn_1:Pkn_4)

df_prior$Type <- rep("Prior", dim(df_prior)[1])
df$Type <- rep("Posterior", dim(df)[1])

df_merged <- rbind(df_prior,df[, c("K","Process_type", "pk", "Type")] )

p_h <- ggplot(df_merged, aes(K, pk, color=Process_type)) +
  geom_bar(aes(linetype=Type), size=0.7, stat="identity", alpha=0.0, position="identity") +
  geom_vline(xintercept=3,  linetype="dashed") + ylab("Density") + xlab(TeX('$K_n$')) +
  theme_minimal() + theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  scale_y_continuous(limits = c(0,1)) +
  scale_color_viridis(discrete="TRUE", begin=0, end=0.9, option="D", name=TeX(sprintf('$\\alpha, \\sigma$')), 
                      labels=unname(TeX(c(sprintf('$\\alpha$ = %.2f, $\\sigma$ = %.1f', alpha_seq[1], sigma_seq[1]),
                                          sprintf('$\\alpha$ = %.2f, $\\sigma$ = %.2f', alpha_seq[1], sigma_seq[2]),
                                          sprintf('$\\alpha$ = %.1f, $\\sigma$ = %.1f', alpha_seq[2], sigma_seq[1]),
                                          sprintf('$\\alpha$ = %.1f, $\\sigma$ = %.2f', alpha_seq[2], sigma_seq[2]))))) +
  scale_linetype_manual(name="Distribution", values=c(1,2), guide=guide_legend(override.aes=list(linetype=c(1, 2), color="black"))) +
  facet_wrap(~Process_type)
plot(p_h)
pdf(paste0(fig_path,"Figure1_PY_thyroid.pdf" ))
plot(p_h)
dev.off()


####################################
######## Plot Univariate PY ########
####################################
final <- loadRData( "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_Slc_PY.RData")
df <- to_df(final)


alpha_seq <- as.numeric(levels(as.factor(df$alpha)))
sigma_seq <- as.numeric(levels(as.factor(df$sigma)))
N <- final[[1]]$N

df <- df%>%group_by(Process_type,alpha,sigma)%>%mutate(pk =pk/sum(pk))
julia <- julia_setup()
julia_library("GibbsTypePriors")

K_ <- max(df$K)
julia_assign("K_bound", K_)
julia_assign("A1", alpha_seq[1])
julia_assign("A2", alpha_seq[2])
julia_assign("S1", sigma_seq[1])
julia_assign("S2", sigma_seq[2])
julia_assign("N", N)

df_prior <- tibble(K = 1:K_, 
                   Pkn_1 = round(julia_eval("Pkn_2PD.(1:K_bound, N, A1, S1)"),3),
                   Pkn_2 = round(julia_eval("Pkn_2PD.(1:K_bound, N, A1, S2)"),3), 
                   Pkn_3 = round(julia_eval("Pkn_2PD.(1:K_bound, N, A2, S1)"),3),
                   Pkn_4 = round(julia_eval("Pkn_2PD.(1:K_bound, N, A2, S2)"),3)) %>% 
  gather(Process_type, pk,Pkn_1:Pkn_4)

df_prior$Type =rep("Prior", dim(df_prior)[1])
df$Type= rep("Posterior", dim(df)[1])

df_merged = rbind(df_prior,df[, c("K","Process_type", "pk", "Type")] )

p_h <- ggplot(df_merged, aes(K, pk, color=Process_type)) +
  geom_bar(aes(linetype=Type),size=0.7, stat="identity", alpha=0.0, position="identity") +
  geom_vline(xintercept=2,  linetype="dashed") + ylab("Density") + xlab(TeX('$K_n$')) +
  theme_minimal() + theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  scale_y_continuous(limits = c(0,1)) +
  scale_color_viridis(discrete="TRUE", begin=0, end=0.9, option="D", name=TeX(sprintf('$\\alpha, \\sigma$')), 
                      labels=unname(TeX(c(sprintf('$\\alpha$ = %.2f, $\\sigma$ = %.1f', alpha_seq[1], sigma_seq[1]),
                                          sprintf('$\\alpha$ = %.2f, $\\sigma$ = %.2f', alpha_seq[1], sigma_seq[2]),
                                          sprintf('$\\alpha$ = %.1f, $\\sigma$ = %.1f', alpha_seq[2], sigma_seq[1]),
                                          sprintf('$\\alpha$ = %.1f, $\\sigma$ = %.2f', alpha_seq[2], sigma_seq[2]))))) +
  scale_linetype_manual(name="Distribution", values=c(1,2), guide=guide_legend(override.aes=list(linetype=c(1, 2), color="black"))) +
  facet_wrap(~Process_type)
plot(p_h)
pdf(paste0(fig_path,"Figure1_PY_Slc.pdf" ))
plot(p_h)
dev.off()