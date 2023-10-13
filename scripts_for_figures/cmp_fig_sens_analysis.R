#rm(list=ls()) 
#setwd("~/Documents/GitHub/BNPconsistency/scripts_for_figures")
## read sources
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Random_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Estimation_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Identification_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Gibbs_sampling_function.R")

require(tidyr)
require(e1071)
require(MCMCpack)
require(mvtnorm)
require(Runuran)
require(flexclust)
library(cowplot)
library(ggplot2)
#---------- B) Specification of the simulation and prior parameters -----------------------------------------------

ds_list<- c("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_20.RData","~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_200.RData",
            "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_2000.RData","~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_20000.RData")


df_sens_1 <-comparison_sens(ds_list, alpha = 0.01,K_ = 10, M_it= 20000 , nburn = 10000, coef_R =1/5)
save(df_sens_1, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens1.RData")

df_sens_2 <-comparison_sens(ds_list, alpha = 0.01,K_ = 10, M_it= 20000 , nburn = 10000, coef_R =1/20)
save(df_sens_2, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens2.RData")

df_sens_3 <-comparison_sens(ds_list, alpha = 1,K_ = 10, M_it= 20000 , nburn = 10000, coef_R =1/5)
save(df_sens_3, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens3.RData")

df_sens_4 <-comparison_sens(ds_list, alpha = 1,K_ = 10, M_it= 20000 , nburn = 10000, coef_R =1/20)
save(df_sens_4, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens4.RData")



#### sensitivity with a fixed Sigma prior
ds_list_short <- c("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_20.RData","~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_200.RData",
                             "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_2000.RData")


df_sens_fixed_1 <-comparison_sens(ds_list_short, alpha = 0.01,K_ = 10, M_it= 20000 , nburn = 10000, coef_R =1, fixS  = TRUE)
save(df_sens_fixed_1, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens_fixed1.RData")

df_sens_fixed_5 <-comparison_sens(ds_list_short, alpha = 0.01,K_ = 10, M_it= 20000 , nburn = 10000, coef_R =1/5, fixS  = TRUE)
save(df_sens_fixed_5, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens_fixed5.RData")


df_sens_fixed_20 <-comparison_sens(ds_list_short, alpha = 0.01,K_ = 10, M_it= 20000 , nburn = 10000, coef_R =1/20, fixS  = TRUE)
save(df_sens_fixed_20, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens_fixed20.RData")


