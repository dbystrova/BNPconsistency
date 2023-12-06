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

df_mult_mtm <- loadRData( "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_thyroid_DMP_MTM.RData")
alpha_l <- as.numeric(levels(as.factor(df_mult_mtm$alpha_level)))

pm <- ggplot(df_mult_mtm, aes(x=c, y=val, color=alpha_level)) + geom_line(aes(linetype=type)) +
  scale_linetype_manual(name="Distribution", values=c("solid", "dotdash"), guide=guide_legend(override.aes=list(color="black"))) +
  scale_color_viridis(discrete="TRUE", begin=0, end=0.9, option="D", name=TeX(sprintf('$\\bar{\\alpha}$')), 
                      labels=unname(TeX(c(sprintf('$\\bar{\\alpha}$ = %.2f', alpha_l[1]),
                                          sprintf('$\\bar{\\alpha}$ = %.1f', alpha_l[2]),
                                          sprintf('$\\bar{\\alpha}$ = %.1f', alpha_l[3]),
                                          sprintf('$\\bar{\\alpha}$ = %1.f', alpha_l[4]))))) +
  theme_minimal() + theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  ylab(expression(tilde(K))) + xlab(TeX('$c$')) + geom_hline(yintercept=3, linetype="dashed", size=0.5, alpha=0.5) +
  facet_wrap(~alpha_level)
plot(pm)
pdf(paste0(fig_path,"Figure_thyroid_DPM_MTM.pdf"))
plot(pm)
dev.off()


################################
######## Univariate DMP ########
################################

df_univ_mtm <- loadRData( "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_Slc_DMP_MTM.RData")
df_univ_mtm <- data.frame(c = df_univ_mtm$с_val, val = df_univ_mtm$val, alpha = df_univ_mtm$alpha, type = df_univ_mtm$type)
alpha_l <- as.numeric(levels(as.factor(df_univ_mtm$alpha)))

pm <- ggplot(df_univ_mtm, aes(x=c, y=val, color=alpha)) + geom_line(aes(linetype=type)) +
  scale_linetype_manual(name="Distribution", values=c("solid", "dotdash"), guide=guide_legend(override.aes=list(color="black"))) +
  scale_color_viridis(discrete="TRUE", begin=0, end=0.9, option="D", name=TeX(sprintf('$\\bar{\\alpha}$')), 
                      labels=unname(TeX(c(sprintf('$\\bar{\\alpha}$ = %.2f', alpha_l[1]),
                                          sprintf('$\\bar{\\alpha}$ = %.1f', alpha_l[2]),
                                          sprintf('$\\bar{\\alpha}$ = %.1f', alpha_l[3]),
                                          sprintf('$\\bar{\\alpha}$ = %1.f', alpha_l[4]))))) +
  theme_minimal() + theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  ylab(expression(tilde(K))) + xlab(TeX('$c$')) + geom_hline(yintercept=2, linetype="dashed", size=0.5, alpha=0.5) +
  facet_wrap(~alpha, scales="free_x")
plot(pm)
pdf(paste0(fig_path,"Figure_Slc_DMP_MTM.pdf"))
plot(pm)
dev.off()