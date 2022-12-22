scripts_for_figures/sim_data/GM_3_20.RData: scripts_for_figures/simdata_generation.R
	Rscript scripts_for_figures/simdata_generation.R
save_for_figures/cmp_fig4.RData: scripts_for_figures/sim_data/GM_3_20.RData scripts_for_figures/cmp_fig4.R
	Rscript scripts_for_figures/cmp_fig4.R
save_for_figures/cmp_fig7.RData: scripts_for_figures/sim_data/GM_3_20.RData scripts_for_figures/cmp_fig7.R
	Rscript scripts_for_figures/cmp_fig7.R
save_for_figures/cmp_fig9.RData: scripts_for_figures/sim_data/GM_3_20.RData scripts_for_figures/cmp_fig9.R
	Rscript scripts_for_figures/cmp_fig9.R
save_for_figures/cmp_fig13.RData: scripts_for_figures/sim_data/GM_3_20.RData scripts_for_figures/cmp_fig13.R
	Rscript scripts_for_figures/cmp_fig13.R
save_for_figures/cmp_fig14.RData: scripts_for_figures/sim_data/GM_3_20.RData scripts_for_figures/cmp_fig14.R
	Rscript scripts_for_figures/cmp_fig14.R
save_for_figures/cmp_fig15.RData: scripts_for_figures/sim_data/GM_3_20.RData scripts_for_figures/cmp_fig15.R
	Rscript scripts_for_figures/cmp_fig15.R
figures/Figure1.pdf: scripts_for_figures/plot_Fig1.R saves_for_figures/cmp_fig4.RData saves_for_figures/cmp_fig7.RData saves_for_figures/cmp_fig13.RData saves_for_figures/cmp_fig14.RData
	Rscript scripts_for_figures/plot_Fig1.R
figures/Figure2.pdf: scripts_for_figures/plt_Figure2.R saves_for_figures/cmp_fig9.RData
	Rscript scripts_for_figures/plot_Fig2.R
figures/Figure3.pdf: scripts_for_figures/plot_Fig3.R saves_for_figures/cmp_fig4.RData saves_for_figures/cmp_fig7.RData saves_for_figures/cmp_fig13.RData saves_for_figures/cmp_fig14.RData
	Rscript scripts_for_figures/plot_Fig3.R
figures/Figure4.pdf: scripts_for_figures/plt_Fig4.R saves_for_figures/cmp_fig15_2.RData
	Rscript scripts_for_figures/plot_Fig4.R
figures: figures/Figure1.pdf \
		figures/Figure2.pdf \
		figures/Figure3.pdf \
		figures/Figure4.pdf
cmp_files: save_for_figures/cmp_fig4.RData \
		   save_for_figures/cmp_fig7.RData \
		   save_for_figures/cmp_fig13.RData \
		   save_for_figures/cmp_fig14.RData 




	
