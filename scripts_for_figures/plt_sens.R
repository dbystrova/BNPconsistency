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
plt_fig1(input_file = "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens1.RData", c_vec =c(0.1, 0.5, 1, 2) , fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure_sens/", Sigma_coef = 1/5 )
# e0 =1
plt_fig1(input_file= "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens2 2.RData", c_vec =c(0.1, 0.5, 1, 2) , fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure_sens/", Sigma_coef = 1/20 )
#e0 = 2.5
plt_fig1(input_file= "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens3.RData", c_vec =c(0.1, 0.5, 1, 2) , fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure_sens/",Sigma_coef = 1/5 )
#e0 = 3
plt_fig1(input_file= "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens4.RData", c_vec =c(0.1, 0.5, 1, 2) , fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure_sens/" ,Sigma_coef = 1/20)


### Plot sensitivity plots with fixed C0

plt_fig1_short(input_file= "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens_fixed1 2.RData", c_vec =c(0.1, 0.5, 1, 2) , fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure_sens/" ,Sigma_coef = 1)


plt_fig1_short(input_file= "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens_fixed5.RData", c_vec =c(0.1, 0.5, 1, 2) , fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure_sens/" ,Sigma_coef = 1/5)



plt_fig1_short(input_file= "~/Documents/GitHub/BNPconsistency/saves_for_figures/cmp_fig_sens_fixed20.RData", c_vec =c(0.1, 0.5, 1, 2) , fig_path= "~/Documents/GitHub/BNPconsistency/figures/Figure_sens/" ,Sigma_coef = 1/20)
