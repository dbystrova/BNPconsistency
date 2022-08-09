
save_for_figures/cmp_fig1.RData: scripts_for_figures/cmp_fig1.R
	Rscript scripts_for_figures/cmp_fig1.R
save_for_figures/cmp_fig2.RData: scripts_for_figures/cmp_fig2.R
	Rscript scripts_for_figures/cmp_fig2.R
save_for_figures/cmp_fig3.RData: scripts_for_figures/cmp_fig3.R
	Rscript scripts_for_figures/cmp_fig3.R
save_for_figures/cmp_fig4.RData: scripts_for_figures/cmp_fig4.R
	Rscript scripts_for_figures/cmp_fig4.R
save_for_figures/cmp_fig5.RData: scripts_for_figures/cmp_fig5.R
	Rscript scripts_for_fgures/cmp_fig5.R	
save_for_figures/cmp_fig6.RData: scripts_for_figures/cmp_fig6.R
	Rscript scripts_for_figures/cmp_fig6.R		
figures/Figure1.pdf: scripts_for_figures/plt_fig1.R saves_for_figures/cmp_fig1.RData 
	Rscript scripts_for_figures/plt_fig1.R
figures/Figure2.pdf: scripts_for_figures/plt_fig2.R saves_for_figures/cmp_fig2.RData 
	Rscript scripts_for_figures/plt_fig2.R
figures/Figure3.pdf: scripts_for_figures/plt_fig3.R saves_for_figures/cmp_fig3.RData 
	Rscript scripts_for_figures/plt_fig3.R
figures/Figure4.pdf: scripts_for_figures/plt_fig4.R saves_for_figures/cmp_fig4.RData 
	Rscript scripts_for_figures/plt_fig4.R
figures/Figure5.pdf: scripts_for_figures/plt_fig5.R saves_for_figures/cmp_fig5.RData 
	Rscript scripts_for_figures/plt_fig5.R
figures/Figure6.pdf: scripts_for_figures/plt_fig6.R saves_for_figures/cmp_fig6.RData 
	Rscript scripts_for_figures/plt_fig6.R	
figures: figures/Figure1.pdf \
		figures/Figure2.pdf \
		figures/Figure3.pdf \
		figures/Figure4.pdf \
		figures/Figure5.pdf \
		figures/Figure6.pdf
	
