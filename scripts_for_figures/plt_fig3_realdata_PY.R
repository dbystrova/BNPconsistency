## read sources
library(tidyverse)
library(latex2exp)
library(viridis)
library(JuliaCall)
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")

fig_path <- "~/Documents/GitHub/BNPconsistency/figures/Figure_real_data/"

#################################
######## Multivariate PY ########
#################################
df_mult_mtm <- loadRData( "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_thyroid_PY_MTM.RData")
names(df_mult_mtm) <- c('c', 'Mean', 'MAP', 'alpha', 'sigma', 'Process_type')
alpha_seq <- as.numeric(levels(as.factor(df_mult_mtm$alpha)))
sigma_seq <- as.numeric(levels(as.factor(df_mult_mtm$sigma)))
df_mult_mtm <- df_mult_mtm %>% gather(key = type, value = val, Mean, MAP) %>%
  group_by(Process_type, alpha, sigma)


pm <- ggplot(df_mult_mtm, aes(x=c, y=val, color=Process_type)) + geom_line(aes(linetype=type)) +
  scale_linetype_manual(name="Distribution", values=c("solid", "dotdash"), guide=guide_legend(override.aes=list(color="black"))) +
  scale_color_viridis(discrete="TRUE", begin=0, end=0.9, option="D", name=TeX(sprintf('$\\alpha, \\sigma$')), 
                      labels=unname(TeX(c(sprintf('$\\alpha$ = %.2f, $\\sigma$ = %.1f', alpha_seq[1], sigma_seq[1]),
                                          sprintf('$\\alpha$ = %.2f, $\\sigma$ = %.2f', alpha_seq[1], sigma_seq[2]),
                                          sprintf('$\\alpha$ = %.1f, $\\sigma$ = %.1f', alpha_seq[2], sigma_seq[1]),
                                          sprintf('$\\alpha$ = %.1f, $\\sigma$ = %.2f', alpha_seq[2], sigma_seq[2])))))+
  theme_minimal() + theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  ylab(expression(tilde(K))) + xlab(TeX('$c$')) + geom_hline(yintercept=3, linetype="dashed", size=0.5, alpha=0.5) +
  facet_wrap(~Process_type)
plot(pm)
pdf(paste0(fig_path,"Figure_thyroid_PY_MTM.pdf"))
plot(pm)
dev.off()


###############################
######## Univariate PY ########
###############################
df_univ_mtm <- loadRData( "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_Slc_PY_MTM.RData")
names(df_univ_mtm) <- c('c', 'Mean', 'MAP', 'alpha', 'sigma', 'Process_type')
alpha_seq <- as.numeric(levels(as.factor(df_univ_mtm$alpha)))
sigma_seq <- as.numeric(levels(as.factor(df_univ_mtm$sigma)))
df_univ_mtm <- df_univ_mtm %>% gather(key = type, value = val, Mean, MAP) %>%
  group_by(Process_type, alpha, sigma)

pm <- ggplot(df_univ_mtm, aes(x=c, y=val, color=Process_type)) + geom_line(aes(linetype=type)) +
  scale_linetype_manual(name="Distribution", values=c("solid", "dotdash"), guide=guide_legend(override.aes=list(color="black"))) +
  scale_color_viridis(discrete="TRUE", begin=0, end=0.9, option="D", name=TeX(sprintf('$\\alpha, \\sigma$')), 
                      labels=unname(TeX(c(sprintf('$\\alpha$ = %.2f, $\\sigma$ = %.1f', alpha_seq[1], sigma_seq[1]),
                                          sprintf('$\\alpha$ = %.2f, $\\sigma$ = %.2f', alpha_seq[1], sigma_seq[2]),
                                          sprintf('$\\alpha$ = %.1f, $\\sigma$ = %.1f', alpha_seq[2], sigma_seq[1]),
                                          sprintf('$\\alpha$ = %.1f, $\\sigma$ = %.2f', alpha_seq[2], sigma_seq[2]))))) +
  theme_minimal() + theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  ylab(expression(tilde(K))) + xlab(TeX('$c$')) + geom_hline(yintercept=2, linetype="dashed", size=0.5, alpha=0.5) +
  facet_wrap(~Process_type)
plot(pm)
pdf(paste0(fig_path,"Figure_Slc_PY_MTM.pdf"))
plot(pm)
dev.off()