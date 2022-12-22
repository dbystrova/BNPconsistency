#rm(list=ls()) 
#setwd("~/Documents/GitHub/BNPconsistency/scripts_for_figures")
## read sources
require(e1071)
require(MASS)
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
source("~/Documents/GitHub/BNPconsistency/scripts_for_figures/plt_Fig4.R")


plt_fig2(input_file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig15.RData", c_v =c(0.1, 0.5, 1, 2),n_list= c(20,200,2000,20000),fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure4_2/",a_g  = 0.1 , b_g =0.1 ) 

plt_fig2(input_file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig15_2.RData", c_v =c(0.1, 0.5, 1, 2),n_list= c(20,200,2000,20000),fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure4_2/",a_g  = 1 , b_g =0.1 ) 



