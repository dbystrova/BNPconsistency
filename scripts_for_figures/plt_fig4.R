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
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Utils_post.R")

#julia_setup(installJulia = FALSE)
#---------- B) Specification of the simulation and prior parameters -----------------------------------------------

fig4 <- loadRData("~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig4.RData")
fig4_ <- fig4$line%>%group_by(Process_type,N)%>%mutate(pkn =density/sum(density))
K_ = max(fig4_$K)
alpha_ = fig4$line$alpha[1]

p = ggplot(fig4_, aes(x=K, colour = fig4_$Process_type)) +
  geom_line(aes(x=K, y = pkn))  +  ylab('')+
  ylab('') + ggtitle(TeX(sprintf('Posterior dist for the number of clusters for $\\N =(%3.f, %2.f, %2.f)$,$\\alpha =%.3f$, \\hat{R} =(%.3f,%.3f,%.3f,,%.3f)',fig4_$N[1],fig4_$N[(max(fig4_$K)+1)],fig4_$N[(2*max(fig4_$K)+1)],fig4_$alpha[1],fig4_$Rh[1],fig4_$Rh[(max(fig4_$K)+1)],fig4_$Rh[(2*max(fig4_$K)+1)],fig4_$Rh[(3*max(fig4_$K)+1)])))+
  theme_minimal() +scale_x_continuous(limits = c(1, max(fig4_$K)), expand = c(0, 0),breaks= c(1,seq(0,max(fig4_$K),length=5)))+
  theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10))+
  scale_color_discrete(name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$N$=%3.f',fig4_$N[1]),sprintf('$N$=%3.f',fig4_$N[(max(fig4_$K)+1)]),sprintf('$N$=%3.f',fig4_$N[(2*max(fig4_$K)+1)]),sprintf('$N$=%3.f',fig4_$N[(3*max(fig4_$K)+1)])))))
p

pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure4_upd.pdf")
plot(p)
dev.off()

julia <- julia_setup()
julia_library("GibbsTypePriors")

julia_assign("K_bound", K_)
julia_assign("N1", 20)
julia_assign("N2", 200)
julia_assign("N3", 2000)
julia_assign("N4", 20000)
# alpha = e0*K
julia_assign("alpha", fig4$line$alpha[1]*K_)
print( fig4$line$alpha[1]*K_)
#a = julia_eval("Pkn_Dirichlet_mult.(1:K_bound,N1, N1, alpha)")

df_prior = tibble(K= 1:K_, 
                  Pkn_1 = round(julia_eval("Pkn_Dirichlet_mult.(1:K_bound,N1, K_bound, alpha)"),3),
                  Pkn_2= round(julia_eval("Pkn_Dirichlet_mult.(1:K_bound,N2, K_bound, alpha)"),3), 
                  Pkn_3 = round(julia_eval("Pkn_Dirichlet_mult.(1:K_bound,N3, K_bound, alpha)"),3),
                  Pkn_4 = round(julia_eval("Pkn_Dirichlet_mult.(1:K_bound,N4, K_bound, alpha)"),3))%>% gather(Process_type, pkn,Pkn_1:Pkn_4)

df_prior$Type =rep("Prior", dim(df_prior)[1])
fig4_$Type= rep("Posterior", dim(fig4_)[1])

df_merged = rbind(df_prior,fig4_[, c("K","Process_type", "pkn", "Type")] )

pkn.labs<- c("Pkn_1","Pkn_2","Pkn_3","Pkn_4")
names(pkn.labs) <- c(paste0("n = ", fig4_$N[1]),paste0("n = ", fig4_$N[max(fig4_$K)+1]),paste0("n = ", fig4_$N[(2*max(fig4_$K)+1)]),paste0("n = ", fig4_$N[(3*max(fig4_$K)+1)]))
pkn_names <- as_labeller(
  c(`Pkn_1` = paste0("n = ", fig4_$N[1]), `Pkn_2` = paste0("n = ", fig4_$N[max(fig4_$K)+1]),`Pkn_3` = paste0("n = ", fig4_$N[(2*max(fig4_$K)+1)]),`Pkn_4` = paste0("n = ", fig4_$N[(3*max(fig4_$K)+1)])))

p <- ggplot(df_merged, aes(K,pkn,color =Process_type))+geom_bar(aes(linetype=Type),size = 0.7, stat="identity",alpha =0.0, position = "identity", fill= "white")+
  geom_vline(xintercept=3,  linetype="dashed")+ylab("Density")+xlab(TeX('$K$'))+
  ggtitle(TeX(sprintf('Posterior distribution for the number of clusters for $\\alpha =%.3f$,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig4_$alpha[1],fig4_$N[1],fig4_$N[(max(fig4_$K)+1)],fig4_$N[(2*max(fig4_$K)+1)],fig4_$N[(3*max(fig4_$K)+1)])))+
  theme_minimal()+ scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$N$=%3.f',fig4_$N[1]),sprintf('$N$=%3.f',fig4_$N[(max(fig4_$K)+1)]),sprintf('$N$=%3.f',fig4_$N[(2*max(fig4_$K)+1)]),sprintf('$N$=%3.f',fig4_$N[(3*max(fig4_$K)+1)])))))+
  scale_x_continuous(breaks = c(1,4,7,10), limits = c(0,11))+
  scale_linetype_manual(name = "Distribution",values = c(1, 2),guide = guide_legend(override.aes = list(linetype = c(1, 2),color = "black") ) )+
  facet_wrap(~Process_type,labeller = pkn_names)
p

pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure4_2_upd.pdf")
plot(p)
dev.off()


p <- ggplot(df_merged, aes(K,pkn,color =Process_type))+geom_bar(aes(linetype=Type),size = 0.7, stat="identity",alpha =0.0, position = "identity", fill= "white")+
  geom_vline(xintercept=3,  linetype="dashed")+ylab("Density")+xlab(TeX('$K_n$'))+
#  ggtitle(TeX(sprintf('Posterior distribution for the number of clusters for $\\alpha =%.3f$,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig4_$alpha[1],fig4_$N[1],fig4_$N[(max(fig4_$K)+1)],fig4_$N[(2*max(fig4_$K)+1)],fig4_$N[(3*max(fig4_$K)+1)])))+
  theme_minimal()+ scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$n$=%3.f',fig4_$N[1]),sprintf('$n$=%3.f',fig4_$N[(max(fig4_$K)+1)]),sprintf('$n$=%3.f',fig4_$N[(2*max(fig4_$K)+1)]),sprintf('$n$=%3.f',fig4_$N[(3*max(fig4_$K)+1)])))))+
  scale_x_continuous(breaks = c(1,4,7,10), limits = c(0,11))+
  scale_linetype_manual(name = "Distribution",values = c(2, 1),breaks=c("Prior","Posterior"),guide = guide_legend(override.aes = list(linetype = c(2, 1),color = "black") ) )+
  facet_wrap(~Process_type,labeller = pkn_names)
p
pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure4_2_no_title_upd.pdf")
plot(p)
dev.off()


E_k<- df_merged %>% group_by(Process_type,Type) %>% summarize(sum(pkn *c(1:K_)))
write.table(E_k, file = "~/Documents/GitHub/BNPconsistency/figures/Prior_Posterior_exp_fig4.csv", sep = ",", col.names = NA,
            qmethod = "double")



fig4_2<- fig4$hist
p <- ggplot(fig4_2, aes(x=density, color =Process_type))+
  geom_histogram(aes(y = ..density..),binwidth = 1, alpha=0.5,position = "identity",fill='white' )+geom_vline(xintercept=3,  linetype="dashed")+
  ggtitle(TeX(sprintf('Posterior distribution for the number of clusters for $\\alpha =%.3f$,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig4_$alpha[1],fig4_$N[1],fig4_$N[(max(fig4_$K)+1)],fig4_$N[(2*max(fig4_$K)+1)],fig4_$N[(3*max(fig4_$K)+1)])))+
  theme_minimal()+scale_x_continuous( limits = c(0,11)) +scale_color_discrete(name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$N$=%3.f',fig4_$N[1]),sprintf('$N$=%3.f',fig4_$N[(max(fig4_$K)+1)]),sprintf('$N$=%3.f',fig4_$N[(2*max(fig4_$K)+1)]),sprintf('$N$=%3.f',fig4_$N[(3*max(fig4_$K)+1)])))))+
  facet_wrap(~Process_type)
p
pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure4_3_upd.pdf")
plot(p)
dev.off()

# Posterior expectation 

## Weights
pkn.labs<- c("W_1","W_2","W_3","W_4")
names(pkn.labs) <- c(paste0("n = ", fig4_$N[1]),paste0("n = ", fig4_$N[max(fig4_$K)+1]),paste0("n = ", fig4_$N[(2*max(fig4_$K)+1)]),paste0("n = ", fig4_$N[(3*max(fig4_$K)+1)]))
pkn_names <- as_labeller(
  c(`W_1` = paste0("n = ", fig4_$N[1]), `W_2` = paste0("n = ", fig4_$N[max(fig4_$K)+1]),`W_3` = paste0("n = ", fig4_$N[(2*max(fig4_$K)+1)]),`W_4` = paste0("n = ", fig4_$N[(3*max(fig4_$K)+1)])))

weights_fig <-fig4$weights 
weights_fig$n <-as.factor(weights_fig$W_val) 

weights_fig_thin<- weights_fig%>%  group_by(K,n, W_,W_val) %>%  filter(row_number() %% 5 == 1)


p <- ggplot(weights_fig_thin, aes(x=K, y=weights, group =K, fill = n )) + ylab("Weights")+xlab(TeX('$K_n$'))+scale_x_continuous(breaks = c(1,4,7,10), limits = c(0,11))+
  geom_boxplot(alpha=0.5) +scale_fill_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D", name = TeX(sprintf('$n$')) ,labels=unname(TeX(c(sprintf('$n$=%3.f',fig4_$N[1]),sprintf('$n$=%3.f',fig4_$N[(max(fig4_$K)+1)]),sprintf('$n$=%3.f',fig4_$N[(2*max(fig4_$K)+1)]),sprintf('$n$=%3.f',fig4_$N[(3*max(fig4_$K)+1)])))))+
  ggtitle(TeX(sprintf('Posterior distribution of the component weights $\\alpha =%.3f$,$\\N =(%2.f,%2.f,%2.f, %2.f) $ ',fig4_$alpha[1],fig4_$N[1],fig4_$N[(max(fig4_$K)+1)],fig4_$N[(2*max(fig4_$K)+1)],fig4_$N[(3*max(fig4_$K)+1)])))+
  theme_minimal()+  facet_wrap(~W_,labeller = pkn_names)
p
pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure4_4.pdf")
plot(p)
dev.off()


c_vec = c(0.1, 0.5, 1, 5)

df_k_20<- df_post(20, c_vec,ind=1, post = fig4  )
df_k_200<- df_post(200, c_vec,ind=2, post = fig4  )
df_k_2000<- df_post(2000, c_vec,ind=3, post = fig4  )
df_k_20000<- df_post(20000, c_vec,ind=4, post = fig4  )

df_fin<- rbind(df_k_20,df_k_200,df_k_2000,df_k_20000)
pkn.labs<- c(20,200,2000,20000)
names(pkn.labs) <- c(paste0("n = ", fig4_$N[1]),paste0("n = ", fig4_$N[max(fig4_$K)+1]),paste0("n = ", fig4_$N[(2*max(fig4_$K)+1)]),paste0("n = ", fig4_$N[(3*max(fig4_$K)+1)]))
pkn_names <- as_labeller(
  c(`20` = paste0("n = ", fig4_$N[1]), `200` = paste0("n = ", fig4_$N[max(fig4_$K)+1]),`2000` = paste0("n = ", fig4_$N[(2*max(fig4_$K)+1)]),`20000` = paste0("n = ", fig4_$N[(3*max(fig4_$K)+1)])))


df_fin_bar = df_fin %>%
  group_by(Process_type,N, pkn) %>%
  summarize(count=n())%>% mutate(pkn_dens =count/sum(count))
K_ = max(fig4_$K)
alpha_ = fig4$line$alpha[1]


p <- ggplot(df_fin_bar, aes(pkn,pkn_dens,color =Process_type, linetype = Process_type))+ geom_bar(aes(linetype=Process_type),size = 0.7, stat="identity",alpha =0.0, position = "identity", fill= "white")+
  geom_vline(xintercept=3,  linetype="dashed")+ylab("Density")+xlab(TeX('$K_n$'))+
  ggtitle(TeX(sprintf('PD for the num of clusters for MTM $\\alpha =%.3f$,$\\c_vec =(%2.1f,%2.1f,%2.1f, %2.1f) $ ',fig4_$alpha[1],c_vec[1],c_vec[2],c_vec[3],c_vec[4])))+
  scale_color_viridis(discrete= "TRUE", begin = 0, end = 0.9,option = "D",name = TeX(sprintf('$c$')) ,labels=unname(TeX(c(sprintf('$c$=%3.1f', c_vec[1]),sprintf('$c$=%3.1f',c_vec[2]),sprintf('$c$=%3.f',c_vec[3]),sprintf('$c$=%3.f',c_vec[4])))))+
  theme_minimal()+ ylab("Density")+xlab(TeX('$K_n$')) +
  scale_linetype_manual(values=c("solid", "dashed","longdash","dotted"),name = TeX(sprintf('$c$')),labels=unname(TeX(c(sprintf('$c$=%3.1f', c_vec[1]),sprintf('$c$=%3.1f',c_vec[2]),sprintf('$c$=%3.f',c_vec[3]),sprintf('$c$=%3.f',c_vec[4]))))) +
  scale_x_continuous(breaks = c(1,4,7,10), limits = c(0,11))+
  facet_wrap(~N)
p

pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure4_5.pdf")
plot(p)
dev.off()

n_vec = c(20,200,2000,20000)

df_MTM_MAP= tibble(c = c_vec,
               Pkn_n1 = as.numeric(colnames(table(df_k_20$Process_type, df_k_20$pkn))[apply(table(df_k_20$Process_type, df_k_20$pkn), 1, which.max)]),
               Pkn_n2=  as.numeric(colnames(table(df_k_200$Process_type, df_k_200$pkn))[apply(table(df_k_200$Process_type, df_k_200$pkn), 1, which.max)]),
               Pkn_n3=  as.numeric(colnames(table(df_k_2000$Process_type, df_k_2000$pkn))[apply(table(df_k_2000$Process_type, df_k_2000$pkn), 1, which.max)]),
               Pkn_n4=  as.numeric(colnames(table(df_k_20000$Process_type, df_k_20000$pkn))[apply(table(df_k_20000$Process_type, df_k_20000$pkn), 1, which.max)]))%>% gather(Process_type, pkn,Pkn_n1:Pkn_n4)


df_MTM_MAP$name<- rep("MAP", dim(df_MTM_MAP)[1])
df_MTM_Mean= tibble(c = c_vec,
                   Pkn_n1 = aggregate( df_k_20$pkn, list(df_k_20$Process_type), FUN=mean)$x,
                   Pkn_n2= aggregate( df_k_200$pkn, list(df_k_20$Process_type), FUN=mean)$x,
                   Pkn_n3=  aggregate( df_k_2000$pkn, list(df_k_20$Process_type), FUN=mean)$x,
                   Pkn_n4=  aggregate( df_k_20000$pkn, list(df_k_20$Process_type), FUN=mean)$x)%>% gather(Process_type, pkn,Pkn_n1:Pkn_n4)


df_MTM_Mean$name<- rep("Mean", dim(df_MTM_MAP)[1])

df_MTM <- rbind(df_MTM_MAP, df_MTM_Mean)

p <- ggplot(df_MTM, aes(c,pkn,color =Process_type)) + geom_point()+
  ggtitle(TeX(sprintf('PD for the num of clusters for MTM $\\alpha =%.3f$,$\\c_vec =(%2.1f,%2.1f,%2.1f, %2.1f) $ ',fig4_$alpha[1],c_vec[1],c_vec[2],c_vec[3],c_vec[4])))+
   theme_minimal()+  facet_wrap(~name)
p

pdf(file="~/Documents/GitHub/BNPconsistency/figures/Figure4_7.pdf")
plot(p)
dev.off()




df_fin_ <- df_fin%>%group_by(Process_type,N)%>%mutate(pkn_dens =pkn_dens/sum(pkn_dens))
K_ = max(fig4_$K)


df = df_fin %>%
  group_by(Process_type,N, pkn) %>%
  summarize(count=n())

