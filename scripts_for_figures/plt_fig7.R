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
library(JuliaCall)
library(viridis)

#---------- B) Specification of the simulation and prior parameters -----------------------------------------------


loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
fig7 <- loadRData("~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig7.RData")
fig7_ <- fig7$line%>%group_by(Process_type,N)%>%mutate(pkn =density/sum(density))

p = ggplot(fig7_, aes(x=K, colour = fig7_$Process_type)) +
  geom_line(aes(x=K, y = pkn))  +  ylab('')+
  ylab('') + ggtitle(TeX(sprintf('Posterior dist for the number of clusters for $\\N =(%3.f, %2.f, %2.f)$,$\\alpha =%2.f$, \\hat{R} =(%.4f,%.4f,%.4f)',fig7_$N[1],fig7_$N[(max(fig7_$K)+1)],fig7_$N[(2*max(fig7_$K)+1)],fig7_$alpha[1],fig7_$Rh[1],fig7_$Rh[(max(fig7_$K)+1)],fig7_$Rh[(2*max(fig7_$K)+1)])))+
  theme_minimal() +scale_x_continuous(limits = c(1, max(fig7_$K)), expand = c(0, 0),breaks= c(1,seq(0,max(fig7_$K),length=5)))+
  theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10))+
  scale_color_discrete(name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$N$=%3.f',fig7_$N[1]),sprintf('$N$=%3.f',fig7_$N[(max(fig7_$K)+1)]),sprintf('$N$=%3.f',fig7_$N[(2*max(fig7_$K)+1)]),sprintf('$N$=%3.f',fig7_$N[(3*max(fig7_$K)+1)])))))
p


pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure7.pdf")
plot(p)
dev.off()


julia <- julia_setup()
julia_library("GibbsTypePriors")


df_prior = tibble(K= 1:10, 
                  Pkn_1 = round(julia_eval("Pkn_Dirichlet_mult.(1:10,20, 10, 1.0)"),3),
                  Pkn_2= round(julia_eval("Pkn_Dirichlet_mult.(1:10,200, 10, 1.0)"),3), 
                  Pkn_3 = round(julia_eval("Pkn_Dirichlet_mult.(1:10,2000, 10, 1.0)"),3),
                  Pkn_4 = round(julia_eval("Pkn_Dirichlet_mult.(1:10,20000, 10, 1.0)"),3))%>% gather(Process_type, pkn,Pkn_1:Pkn_4)

df_prior$Type =rep("Prior", dim(df_prior)[1])
fig7_$Type= rep("Posterior", dim(fig7_)[1])

df_merged = rbind(df_prior,fig7_[, c("K","Process_type", "pkn", "Type")] )

pkn.labs<- c("Pkn_1","Pkn_2","Pkn_3","Pkn_4")
names(pkn.labs) <- c(paste0("N = ", fig7_$N[1]),paste0("N = ", fig7_$N[max(fig7_$K)+1]),paste0("N = ", fig7_$N[(2*max(fig7_$K)+1)]),paste0("N = ", fig7_$N[(3*max(fig7_$K)+1)]))
pkn_names <- as_labeller(
  c(`Pkn_1` = paste0("N = ", fig7_$N[1]), `Pkn_2` = paste0("N = ", fig7_$N[max(fig7_$K)+1]),`Pkn_3` = paste0("N = ", fig7_$N[(2*max(fig7_$K)+1)]),`Pkn_4` = paste0("N = ", fig7_$N[(3*max(fig7_$K)+1)])))

p <- ggplot(df_merged, aes(K,pkn,color =Process_type))+geom_bar(aes(linetype=Type),size = 0.7, stat="identity",alpha =1, position = "identity", fill= "white")+
  geom_vline(xintercept=3,  linetype="dashed")+ylab("Density")+xlab(TeX('$K$'))+
  ggtitle(TeX(sprintf('Posterior distribution for the number of clusters for $\\alpha =%.3f$,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig7_$alpha[1],fig7_$N[1],fig7_$N[(max(fig7_$K)+1)],fig7_$N[(2*max(fig7_$K)+1)],fig7_$N[(3*max(fig7_$K)+1)])))+
  theme_minimal()+ scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$N$=%3.f',fig7_$N[1]),sprintf('$N$=%3.f',fig7_$N[(max(fig7_$K)+1)]),sprintf('$N$=%3.f',fig7_$N[(2*max(fig7_$K)+1)]),sprintf('$N$=%3.f',fig7_$N[(3*max(fig7_$K)+1)])))))+
  scale_x_continuous(breaks = c(1,3,6,9), limits = c(0,11))+
  scale_linetype_manual(name = "Distribution",values = c(1, 2),guide = guide_legend(override.aes = list(linetype = c(1, 2),color = "black") ) )+
  facet_wrap(~Process_type,labeller = pkn_names)
p

pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure7_2.pdf")
plot(p)
dev.off()


fig7_2 <- fig7$hist
p <- ggplot(fig7_2, aes(x=density, color =Process_type))+
  geom_histogram(binwidth = 1, alpha=0.5,position = "identity",fill='white' )+geom_vline(xintercept=3,  linetype="dashed")+
  ggtitle(TeX(sprintf('Posterior distribution for the number of clusters for $\\alpha =%.3f$,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig7_$alpha[1],fig7_$N[1],fig7_2$N[(max(fig7_$K)+1)],fig7_$N[(2*max(fig7_$K)+1)],fig7_$N[(3*max(fig7_$K)+1)])))+
  theme_minimal()+scale_color_discrete(name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$N$=%3.f',fig7_$N[1]),sprintf('$N$=%3.f',fig7_$N[(max(fig7_$K)+1)]),sprintf('$N$=%3.f',fig7_$N[(2*max(fig7_$K)+1)]),sprintf('$N$=%3.f',fig7_$N[(3*max(fig7_$K)+1)])))))+
  facet_wrap(~Process_type)
p
pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure7_3.pdf")
plot(p)
dev.off()

