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
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/MTM.R")

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
#---------- B) Specification of the simulation and prior parameters -----------------------------------------------

fig5 <- loadRData("~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig5.RData")
fig5_ <- fig5$line%>%group_by(Process_type,N)%>%mutate(pkn =density/sum(density))

K_ = max(fig5_$K)
#E_k<- fig5_ %>% group_by(Process_type) %>% summarize(sum(pkn *c(1:K_)))
#write.table(E_k, file = "Posterior_exp_fig5.csv", sep = ",", col.names = NA,
#            qmethod = "double")



p = ggplot(fig5_, aes(x=K, colour = fig5_$Process_type)) +
  geom_line(aes(x=K, y = pkn))  +  ylab('')+
  ylab('') + ggtitle(TeX(sprintf('Posterior dist for the number of clusters for $\\N =(%3.f, %2.f, %2.f)$,$\\alpha =%.2f$, \\hat{R} =(%.3f,%.3f,%.3f,%.3f)',fig5_$N[1],fig5_$N[(max(fig5_$K)+1)],fig5_$N[(2*max(fig5_$K)+1)],fig5_$alpha[1],fig5_$Rh[1],fig5_$Rh[(max(fig5_$K)+1)],fig5_$Rh[(2*max(fig5_$K)+1)],fig5_$Rh[(2*max(fig5_$K)+1)])))+
  theme_minimal() +scale_x_continuous(limits = c(1, max(fig5_$K)), expand = c(0, 0),breaks= c(1,seq(0,max(fig5_$K),length=5)))+
  theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10))+geom_vline(xintercept=3,  linetype="dashed")+
  scale_color_discrete(name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$N$=%3.f',fig5_$N[1]),sprintf('$N$=%3.f',fig5_$N[(max(fig5_$K)+1)]),sprintf('$N$=%3.f',fig5_$N[(2*max(fig5_$K)+1)]),sprintf('$N$=%3.f',fig5_$N[(3*max(fig5_$K)+1)])))))
p


pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure5_upd.pdf")
plot(p)
dev.off()

julia <- julia_setup()
julia_library("GibbsTypePriors")


df_prior = tibble(K= 1:10, 
                  Pkn_1 = round(julia_eval("Pkn_Dirichlet_mult.(1:10,20, 10, 5.0)"),3),
                  Pkn_2= round(julia_eval("Pkn_Dirichlet_mult.(1:10,200, 10, 5.0)"),3), 
                  Pkn_3 = round(julia_eval("Pkn_Dirichlet_mult.(1:10,2000, 10, 5.0)"),3),
                  Pkn_4 = round(julia_eval("Pkn_Dirichlet_mult.(1:10,20000, 10, 5.0)"),3))%>% gather(Process_type, pkn,Pkn_1:Pkn_4)

df_prior$Type =rep("Prior", dim(df_prior)[1])
fig5_$Type= rep("Posterior", dim(fig5_)[1])

df_merged = rbind(df_prior,fig5_[, c("K","Process_type", "pkn", "Type")] )

pkn.labs<- c("Pkn_1","Pkn_2","Pkn_3","Pkn_4")
names(pkn.labs) <- c(paste0("N = ", fig5_$N[1]),paste0("N = ", fig5_$N[max(fig5_$K)+1]),paste0("N = ", fig5_$N[(2*max(fig5_$K)+1)]),paste0("N = ", fig5_$N[(3*max(fig5_$K)+1)]))
pkn_names <- as_labeller(
  c(`Pkn_1` = paste0("N = ", fig5_$N[1]), `Pkn_2` = paste0("N = ", fig5_$N[max(fig5_$K)+1]),`Pkn_3` = paste0("N = ", fig5_$N[(2*max(fig5_$K)+1)]),`Pkn_4` = paste0("N = ", fig5_$N[(3*max(fig5_$K)+1)])))

p <- ggplot(df_merged, aes(K,pkn,color =Process_type))+geom_bar(aes(linetype=Type),size = 0.7, stat="identity",alpha =0.0, position = "identity")+
  geom_vline(xintercept=3,  linetype="dashed")+ylab("Density")+xlab(TeX('$K$'))+
  ggtitle(TeX(sprintf('Posterior distribution for the number of clusters for $\\alpha =%.3f$,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig5_$alpha[1],fig5_$N[1],fig5_$N[(max(fig5_$K)+1)],fig5_$N[(2*max(fig5_$K)+1)],fig5_$N[(3*max(fig5_$K)+1)])))+
  theme_minimal()+ scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$N$=%3.f',fig5_$N[1]),sprintf('$N$=%3.f',fig5_$N[(max(fig5_$K)+1)]),sprintf('$N$=%3.f',fig5_$N[(2*max(fig5_$K)+1)]),sprintf('$N$=%3.f',fig5_$N[(3*max(fig5_$K)+1)])))))+
  scale_x_continuous(breaks = c(1,3,6,9), limits = c(0,11))+
  scale_linetype_manual(name = "Distribution",values = c(1, 2),guide = guide_legend(override.aes = list(linetype = c(1, 2),color = "black") ) )+
  facet_wrap(~Process_type,labeller = pkn_names)
p
pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure5_2_upd.pdf")
plot(p)
dev.off()

E_k<- df_merged %>% group_by(Process_type,Type) %>% summarize(sum(pkn *c(1:K_)))
write.table(E_k, file = "~/Documents/GitHub/BNPconsistency/figures/Prior_Posterior_exp_fig5.csv", sep = ",", col.names = NA,
            qmethod = "double")


fig5_2<- fig5$hist

fig5_2 <- fig5$hist
p <- ggplot(fig5_2, aes(x=density, color =Process_type))+
  geom_histogram(binwidth = 1, alpha=0.2,position = "identity", fill = "white")+geom_vline(xintercept=3,  linetype="dashed")+
  ylab('') + ggtitle(TeX(sprintf('PD for the number of clusters for $\\N =(%3.f, %2.f, %2.f)$,$\\alpha =%.2f$, \\hat{R} =(%.4f,%.4f,%.4f)',fig5_$N[1],fig5_$N[(max(fig5_$K)+1)],fig5_$N[(2*max(fig5_$K)+1)],fig5_$alpha[1],fig5_$Rh[1],fig5_$Rh[(max(fig5_$K)+1)],fig5_$Rh[(2*max(fig5_$K)+1)])))+
  theme_minimal()+scale_color_discrete(name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$N$=%3.f',fig5_$N[1]),sprintf('$N$=%3.f',fig5_$N[(max(fig5_$K)+1)]),sprintf('$N$=%3.f',fig5_$N[(2*max(fig5_$K)+1)]),sprintf('$N$=%3.f',fig5_$N[(3*max(fig5_$K)+1)])))))+
  facet_wrap(~Process_type)
p
pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure5_3_upd.pdf")
plot(p)
dev.off()



weights_fig <-fig5$weights 
weights_fig$n <-as.factor(weights_fig$W_val) 
p <- ggplot(weights_fig, aes(x=K, y=weights, group =K, fill = n )) + 
  geom_boxplot(alpha=0.5) +scale_fill_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$n$=%3.f',fig5_$N[1]),sprintf('$n$=%3.f',fig5_$N[(max(fig5_$K)+1)]),sprintf('$n$=%3.f',fig5_$N[(2*max(fig5_$K)+1)]),sprintf('$n$=%3.f',fig5_$N[(3*max(fig5_$K)+1)])))))+
  ggtitle(TeX(sprintf('Posterior distribution of the component weights $\\alpha =%.3f$,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig5_$alpha[1],fig5_$N[1],fig5_$N[(max(fig5_$K)+1)],fig5_$N[(2*max(fig5_$K)+1)],fig5_$N[(3*max(fig5_$K)+1)])))+
  theme_minimal()+  facet_wrap(~W_val)
p
pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure5_4.pdf")
plot(p)
dev.off()


df_post<- function(n, c_vec,ind, post){
  df_post_k= tibble(it= 1:dim(post$eta[[ind]])[1],
                    Pkn_c1 = round(apply_MTM(post$eta[[ind]],post$mu[[ind]],post$sig[[ind]], c =c_vec[1], n=n),3),
                    Pkn_c2= round(apply_MTM(post$eta[[ind]],post$mu[[ind]],post$sig[[ind]], c = c_vec[2], n=n),3), 
                    Pkn_c3 = round(apply_MTM(post$eta[[ind]],post$mu[[ind]],post$sig[[ind]], c =c_vec[3], n=n),3),
                    Pkn_c4 = round(apply_MTM(post$eta[[ind]],post$mu[[ind]],post$sig[[ind]], c = c_vec[4], n=n),3))%>% gather(Process_type, pkn,Pkn_c1:Pkn_c4)
  
  df_post_k$N = rep(n,dim(df_post_k)[1])
  return(df_post_k)
}

c_vec = c(0.1, 0.5, 1, 5)

df_k_20<- df_post(20, c_vec=c_vec,ind=1, post = fig5  )
df_k_200<- df_post(200, c_vec = c_vec,ind=2, post = fig5  )
df_k_2000<- df_post(2000, c_vec=c_vec,ind=3, post = fig5  )
df_k_20000<- df_post(20000,c_vec=c_vec,ind=4, post = fig5  )
df_fin<- rbind(df_k_20,df_k_200,df_k_2000,df_k_20000)



p <- ggplot(df_fin, aes(pkn,color =Process_type)) + geom_histogram(aes(y = ..density..),binwidth = 1, alpha=0.5,position = "identity",fill='white' )+geom_vline(xintercept=3,  linetype="dashed")+
  ggtitle(TeX(sprintf('PD for the num of clusters for MTM $\\alpha =%.3f$,$\\c_vec =(%2.1f,%2.1f,%2.1f, %2.1f) $ ',fig5_$alpha[1],c_vec[1],c_vec[2],c_vec[3],c_vec[4])))+
  scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D",name = TeX(sprintf('$c$')) ,labels=unname(TeX(c(sprintf('$c$=%3.1f', c_vec[1]),sprintf('$c$=%3.1f',c_vec[2]),sprintf('$c$=%3.f',c_vec[3]),sprintf('$c$=%3.f',c_vec[4])))))+
  theme_minimal()+ 
  scale_x_continuous(breaks = c(1,3,6,9), limits = c(0,11))+
  facet_wrap(~N)
p

pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure5_6.pdf")
plot(p)
dev.off()
