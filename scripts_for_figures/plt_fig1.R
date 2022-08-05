#rm(list=ls()) 
#setwd("~/Documents/GitHub/BNPconsistency/scripts_for_figures")
## read sources
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
fig1 <- loadRData("../saves_for_figures/cmp_fig1.RData")
fig1_ <- fig1%>%group_by(Process_type,alpha)%>%mutate(pkn =density/sum(density))

p = ggplot(fig1_, aes(x=K, y = pkn, colour = fig1_$Process_type)) +
  geom_line(aes(x=K, y = pkn))  +
  ylab('') + ggtitle(TeX(sprintf('Posterior distribution for the number of clusters for $\\alpha =(%.3f, %.2f, %.2f)$,$\\N =%2.f$ ',fig1$alpha[1],fig1$alpha[(max(fig1$K)+1)],fig1$alpha[(2*max(fig1$K)+1)],fig1$N[1])))+
  theme_minimal() +scale_x_continuous(limits = c(1, max(fig1$K)), expand = c(0, 0),breaks= c(1,seq(0,max(fig1$K),length=5)))+
  theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10))+
  scale_color_discrete(name = TeX(sprintf('$\\alpha$')) ,labels=unname(TeX(c(sprintf('$\\alpha$=%.3f',fig1$alpha[1]),sprintf('$\\alpha$=%.3f',fig1$alpha[(max(fig1$K)+1)]),sprintf('$\\alpha$=%.3f',fig1$alpha[(2*max(fig1$K)+1)])))))
p
                     

pdf(file="../figures/Figure1.pdf")
plot(p)
dev.off()
