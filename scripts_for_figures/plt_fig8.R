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
fig8 <- loadRData("~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig8.RData")
fig8_ <- fig8$line%>%group_by(Process_type,N)%>%mutate(pkn =density/sum(density))

p = ggplot(fig8_, aes(x=K, y = pkn, colour = fig8_$Process_type)) +
  geom_line(aes(x=K, y = pkn))  +ylab('') + 
   ggtitle(TeX(sprintf('PD for the number of clusters for $\\N =(%3.f, %2.f, %2.f, %2.f)$,$\\alpha =%.2f$, \\hat{R} =(%.4f,%.4f,%.4f)',fig8_$N[1],fig8_$N[(max(fig8_$K)+1)],fig8_$N[(2*max(fig8_$K)+1)],fig8_$N[(3*max(fig8_$K)+1)],fig8_$alpha[1],fig8_$Rh[1],fig8_$Rh[(max(fig8_$K)+1)],fig8_$Rh[(2*max(fig8_$K)+1)])))+
  theme_minimal() +scale_x_continuous(limits = c(1, max(fig8_$K)), expand = c(0, 0),breaks= c(1,seq(0,max(fig8_$K),length=5)))+
  theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10))+
  scale_color_discrete(name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$N$=%3.f',fig8_$N[1]),sprintf('$N$=%3.f',fig8_$N[(max(fig8_$K)+1)]),sprintf('$N$=%3.f',fig8_$N[(2*max(fig8_$K)+1)]),sprintf('$N$=%3.f',fig8_$N[(3*max(fig8_$K)+1)])))))
p

pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure8.pdf")
plot(p)
dev.off()

fig8_2 = fig8$hist
p <- ggplot(fig8_2, aes(x=density, color =Process_type))+
  geom_histogram(binwidth = 1, alpha=0.5,position = "identity",fill='white' )+geom_vline(xintercept=3,  linetype="dashed")+
  ggtitle(TeX(sprintf('PD for the number of clusters for $\\N =(%3.f, %2.f, %2.f, %2.f)$,$\\alpha =%.2f$, \\hat{R} =(%.4f,%.4f,%.4f)',fig8_2$N[1],fig8_2$N[(max(fig8_2$K)+1)],fig8_2$N[(2*max(fig8_2$K)+1)],fig8_2$N[(3*max(fig8_2$K)+1)],fig8_2$alpha[1],fig8_2$Rh[1],fig8_2$Rh[(max(fig8_2$K)+1)],fig8_2$Rh[(2*max(fig8_2$K)+1)])))+
  theme_minimal()+scale_color_discrete(name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$N$=%3.f',fig8_2$N[1]),sprintf('$N$=%3.f',fig8_2$N[(max(fig8_2$K)+1)]),sprintf('$N$=%3.f',fig8_2$N[(2*max(fig8_2$K)+1)]),sprintf('$N$=%3.f',fig8_2$N[(3*max(fig8_2$K)+1)])))))+
  facet_wrap(~Process_type)
p
pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure8_2.pdf")
plot(p)
dev.off()

