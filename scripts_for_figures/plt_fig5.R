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
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
#---------- B) Specification of the simulation and prior parameters -----------------------------------------------

fig5 <- loadRData("~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig5.RData")
fig5_ <- fig5$line%>%group_by(Process_type,N)%>%mutate(pkn =density/sum(density))


p = ggplot(fig5_, aes(x=K, colour = fig5_$Process_type)) +
  geom_line(aes(x=K, y = pkn))  +  ylab('')+
  ylab('') + ggtitle(TeX(sprintf('Posterior dist for the number of clusters for $\\N =(%3.f, %2.f, %2.f)$,$\\alpha =%.2f$, \\hat{R} =(%.3f,%.3f,%.3f,%.3f)',fig5_$N[1],fig5_$N[(max(fig5_$K)+1)],fig5_$N[(2*max(fig5_$K)+1)],fig5_$alpha[1],fig5_$Rh[1],fig5_$Rh[(max(fig5_$K)+1)],fig5_$Rh[(2*max(fig5_$K)+1)],fig5_$Rh[(2*max(fig5_$K)+1)])))+
  theme_minimal() +scale_x_continuous(limits = c(1, max(fig5_$K)), expand = c(0, 0),breaks= c(1,seq(0,max(fig5_$K),length=5)))+
  theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10))+geom_vline(xintercept=3,  linetype="dashed")+
  scale_color_discrete(name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$N$=%3.f',fig5_$N[1]),sprintf('$N$=%3.f',fig5_$N[(max(fig5_$K)+1)]),sprintf('$N$=%3.f',fig5_$N[(2*max(fig5_$K)+1)]),sprintf('$N$=%3.f',fig5_$N[(3*max(fig5_$K)+1)])))))
p


pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure5.pdf")
plot(p)
dev.off()


fig5_2 <- fig5$hist
p <- ggplot(fig5_2, aes(x=density, color =Process_type))+
  geom_histogram(binwidth = 1, alpha=0.2,position = "identity", fill = "white")+geom_vline(xintercept=3,  linetype="dashed")+
  ylab('') + ggtitle(TeX(sprintf('PD for the number of clusters for $\\N =(%3.f, %2.f, %2.f)$,$\\alpha =%.2f$, \\hat{R} =(%.4f,%.4f,%.4f)',fig5_$N[1],fig5_$N[(max(fig5_$K)+1)],fig5_$N[(2*max(fig5_$K)+1)],fig5_$alpha[1],fig5_$Rh[1],fig5_$Rh[(max(fig5_$K)+1)],fig5_$Rh[(2*max(fig5_$K)+1)])))+
  theme_minimal()+scale_color_discrete(name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$N$=%3.f',fig5_$N[1]),sprintf('$N$=%3.f',fig5_$N[(max(fig5_$K)+1)]),sprintf('$N$=%3.f',fig5_$N[(2*max(fig5_$K)+1)]),sprintf('$N$=%3.f',fig5_$N[(3*max(fig5_$K)+1)])))))+
  facet_wrap(~Process_type)
p
pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure5_2.pdf")
plot(p)
dev.off()
