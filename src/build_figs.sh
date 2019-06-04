# makes complicated SVG figs that rely on access to the EiaD sqlite file
# have to hand convert SVG to png for manuscript

cd ~/git/eyeIntegration_v1_app_manuscript
Rscript src/data_pull.R
Rscript src/retina_developmental_heatmaps.R
Rscript src/tsne_and_pca_calcs.R
Rscript figures_and_tables/zoom_heatmap_figure.R