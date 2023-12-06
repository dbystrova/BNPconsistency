## read sources
library(tidyverse)
library(latex2exp)
library(viridis)
library(JuliaCall)
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")

fig_path <- "~/Documents/GitHub/BNPconsistency/figures/Figure_real_data/"
##################################
######## Multivariate DMP ########
##################################

final_mult <- loadRData( "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_thyroid_DMP.RData")
alpha_l <- as.numeric(levels(as.factor(final_mult$line$Al)))
N <- final_mult$hist$N[1]


fig_df_mut <- final_mult$line%>%group_by(Process_type,Al)%>%mutate(pkn =density/sum(density))
julia <- julia_setup()
julia_library("GibbsTypePriors")

K_ <- 10
julia_assign("K_bound", K_)
julia_assign("A1", alpha_l[1]*K_)
julia_assign("A2", alpha_l[2]*K_)
julia_assign("A3", alpha_l[3]*K_)
julia_assign("A4", alpha_l[4]*K_)
julia_assign("N", N)

df_prior <- tibble(K= 1:K_, 
                  Pkn_1 = round(julia_eval("Pkn_Dirichlet_mult.(1:10, N, 10, A1)"),3),
                  Pkn_2 = round(julia_eval("Pkn_Dirichlet_mult.(1:10, N, 10, A2)"),3), 
                  Pkn_3 = round(julia_eval("Pkn_Dirichlet_mult.(1:10, N, 10, A3)"),3),
                  Pkn_4 = round(julia_eval("Pkn_Dirichlet_mult.(1:10, N, 10, A4)"),3))%>% gather(Process_type, pkn,Pkn_1:Pkn_4)

df_prior$Type <- rep("Prior", dim(df_prior)[1])
fig_df_mut$Type <- rep("Posterior", dim(fig_df_mut)[1])

df_merged <- rbind(df_prior,fig_df_mut[, c("K","Process_type", "pkn", "Type")] )


p_h <- ggplot(df_merged, aes(K, pkn, color=Process_type)) +
  geom_bar(aes(linetype=Type),size=0.7, stat="identity", alpha=0.0, position="identity") +
  geom_vline(xintercept=3,  linetype="dashed") + ylab("Density") + xlab(TeX('$K_n$')) +
  theme_minimal() + theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  scale_color_viridis(discrete="TRUE", begin=0, end=0.9, option="D", name=TeX(sprintf('$\\bar{\\alpha}$')), 
                      labels=unname(TeX(c(sprintf('$\\bar{\\alpha}$ = %.2f', alpha_l[1]),
                                          sprintf('$\\bar{\\alpha}$ = %.1f', alpha_l[2]),
                                          sprintf('$\\bar{\\alpha}$ = %.1f', alpha_l[3]),
                                          sprintf('$\\bar{\\alpha}$ = %1.f', alpha_l[4]))))) +
  scale_x_continuous(breaks=c(1,3,7,10), limits=c(0,11)) +
  scale_linetype_manual(name="Distribution", values=c(1,2), guide=guide_legend(override.aes=list(linetype=c(1, 2), color="black"))) +
  facet_wrap(~Process_type)
p_h
pdf(paste0(fig_path,"Figure1_DMP_thyroid.pdf" ))
plot(p_h)
dev.off()




################################
######## Univariate DMP ########
################################

final_univ <- loadRData( "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_Slc_DMP.RData")
alpha_l <- as.numeric(levels(as.factor(final_univ$line$Al)))
N <- final_univ$hist$N[1]

fig_df_mut <- final_univ$line%>%group_by(Process_type,Al)%>%mutate(pkn =density/sum(density))
julia <- julia_setup()
julia_library("GibbsTypePriors")

K_ <- 10
julia_assign("K_bound", K_)
julia_assign("A1", alpha_l[1]*K_)
julia_assign("A2", alpha_l[2]*K_)
julia_assign("A3", alpha_l[3]*K_)
julia_assign("A4", alpha_l[4]*K_)
julia_assign("N", N)

df_prior <- tibble(K= 1:K_, 
                  Pkn_1 = round(julia_eval("Pkn_Dirichlet_mult.(1:10,215, 10, A1)"),3),
                  Pkn_2= round(julia_eval("Pkn_Dirichlet_mult.(1:10,215, 10, A2)"),3), 
                  Pkn_3 = round(julia_eval("Pkn_Dirichlet_mult.(1:10,215, 10, A3)"),3),
                  Pkn_4 = round(julia_eval("Pkn_Dirichlet_mult.(1:10,215, 10, A4)"),3))%>% gather(Process_type, pkn,Pkn_1:Pkn_4)

df_prior$Type <- rep("Prior", dim(df_prior)[1])
fig_df_mut$Type <- rep("Posterior", dim(fig_df_mut)[1])

df_merged <- rbind(df_prior,fig_df_mut[, c("K","Process_type", "pkn", "Type")] )

p_h <- ggplot(df_merged, aes(K, pkn, color=Process_type)) +
  geom_bar(aes(linetype=Type),size=0.7, stat="identity", alpha=0.0, position="identity") +
  geom_vline(xintercept=2,  linetype="dashed") + ylab("Density") + xlab(TeX('$K_n$')) +
  theme_minimal() + theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  scale_color_viridis(discrete="TRUE", begin=0, end=0.9, option="D", name=TeX(sprintf('$\\bar{\\alpha}$')), 
                      labels=unname(TeX(c(sprintf('$\\bar{\\alpha}$ = %.2f', alpha_l[1]),
                                          sprintf('$\\bar{\\alpha}$ = %.1f', alpha_l[2]),
                                          sprintf('$\\bar{\\alpha}$ = %.1f', alpha_l[3]),
                                          sprintf('$\\bar{\\alpha}$ = %1.f', alpha_l[4]))))) +
  scale_x_continuous(breaks=c(1,2,4,6,8,10), limits=c(0,11)) +
  scale_linetype_manual(name="Distribution", values=c(1,2), guide=guide_legend(override.aes=list(linetype=c(1, 2), color="black"))) +
  facet_wrap(~Process_type)
p_h
pdf(paste0(fig_path,"Figure1_Slc_DMP.pdf" ))
plot(p_h)
dev.off()