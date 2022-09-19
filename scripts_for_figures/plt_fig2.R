#rm(list=ls()) 
require(e1071)
require(mclust)
require(MASS)
require(bayesm)
require(MCMCpack)
require(mvtnorm)
require(Runuran)
require(flexclust)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(latex2exp)
require(tidyr)
library(dplyr)

#---------- B) Specification of the simulation and prior parameters -----------------------------------------------


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
fig2 <- loadRData("~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig2.RData")
fig2_ <- fig2%>%group_by(Process_type,alpha)%>%mutate(pkn =density/sum(density))

p = ggplot(fig2_, aes(x=K, y = pkn, colour = fig2_$Process_type)) +
  geom_line(aes(x=K, y = pkn))  +
  ylab('') + ggtitle(TeX(sprintf('Posterior dist for the number of clusters for $\\alpha =(%.3f, %.2f, %.2f)$,$\\N =%2.f$, \\hat{R} =(%.4f,%.4f,%.4f)',fig2_$alpha[1],fig2_$alpha[(max(fig2_$K)+1)],fig2_$alpha[(2*max(fig2_$K)+1)],fig2_$N[1],fig2_$Rh[1],fig2_$Rh[(max(fig2_$K)+1)],fig2_$Rh[(2*max(fig2_$K)+1)])))+
  theme_minimal() +scale_x_continuous(limits = c(1, max(fig2_$K)), expand = c(0, 0),breaks= c(1,seq(0,max(fig2_$K),length=5)))+
  theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10))+
  scale_color_discrete(name = TeX(sprintf('$\\alpha$')) ,labels=unname(TeX(c(sprintf('$\\alpha$=%.3f',fig2_$alpha[1]),sprintf('$\\alpha$=%.3f',fig2_$alpha[(max(fig2_$K)+1)]),sprintf('$\\alpha$=%.3f',fig2_$alpha[(2*max(fig2_$K)+1)])))))
p


pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure2.pdf")
plot(p)
dev.off()
