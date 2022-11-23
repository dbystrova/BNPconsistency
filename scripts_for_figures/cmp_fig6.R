rm(list=ls()) 
#setwd("~/Documents/GitHub/BNPconsistency/scripts_for_figures")
## read sources
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Random_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Estimation_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Code_SP_Mix/Identification_SpMix.R")
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/Gibbs_sampling_function.R")

require(tidyr)
require(e1071)
require(MASS)
require(MCMCpack)
require(mvtnorm)
require(Runuran)
require(flexclust)
library(cowplot)
library(ggplot2)
#---------- B) Specification of the simulation and prior parameters -----------------------------------------------


ds_list<- c("~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_20.RData","~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_200.RData",
            "~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_2000.RData","~/Documents/GitHub/BNPconsistency/scripts_for_figures/sim_data/GM_3_20000.RData")


df_6 <-comparison_n(ds_list, alpha = 0.9,K_ = 10, M_it= 10000 , nburn = 2000)
save(df_6, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig6_test.RData")

#save(df_6, file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig6.RData")


