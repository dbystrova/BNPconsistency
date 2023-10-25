rm(list=ls()) 
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
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/plt_Figure1.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/plt_Fig1_short.R")

# e0 = 0.01
plt_fig1(input_file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens1.RData", c_vec =c(0.1, 0.5, 1, 2) , fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure_sensitivity/", Sigma_coef = 1/5 , MTM = FALSE)
# e0 =1
plt_fig1(input_file= "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens2 2.RData", c_vec =c(0.1, 0.5, 1, 2) , fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure_sensitivity/", Sigma_coef = 1/20, MTM = FALSE )
#e0 = 2.5
plt_fig1(input_file= "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens3.RData", c_vec =c(0.1, 0.5, 1, 2) , fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure_sensitivity/",Sigma_coef = 1/5 , MTM = FALSE)
#e0 = 3
plt_fig1(input_file= "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens4.RData", c_vec =c(0.1, 0.5, 1, 2) , fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure_sensitivity/" ,Sigma_coef = 1/20, MTM = FALSE)


#plt_fig1_short(input_file= "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens_fixed1 2.RData", c_vec =c(0.1, 0.5, 1, 2) , fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure_sens/" ,Sigma_coef = 1)
#plt_fig1_short(input_file= "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens_fixed5 2.RData", c_vec =c(0.1, 0.5, 1, 2) , fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure_sens/" ,Sigma_coef = 1/5)
#plt_fig1_short(input_file= "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens_fixed20.RData", c_vec =c(0.1, 0.5, 1, 2) , fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure_sens/" ,Sigma_coef = 1/20)
#plt_fig1_short(input_file= "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens_fixed05.RData", c_vec =c(0.1, 0.5, 1, 2) , fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure_sens/" ,Sigma_coef = 5)


#Convergence results 
df_1 =  loadRData("~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens1.RData")
df_2 =  loadRData("~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens2 2.RData")
df_3 =  loadRData("~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens3.RData")
df_4 =  loadRData("~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens4.RData")

df_fix  =  loadRData("~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens_ds20.RData")

Rh2 = df_2$line$Rh
Rh2[1:10] = df_fix$line$Rh
Rh_table = tibble(K= 1:40,Rh_1 = df_1$line$Rh, Rh_2 = Rh2, Rh_3 = df_3$line$Rh , Rh_4 = df_4$line$Rh)

write.table(Rh_table, file =paste0(fig_path,"Convergence_sensitivity_analysis.csv"), sep = ",", col.names = NA,qmethod = "double")

