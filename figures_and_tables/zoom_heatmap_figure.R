load('data/data_pull.Rdata')
library(tidyverse)
library(ggforce)
library(dbscan)
library(cowplot)

retina_plus_group <- tsne_50 %>% filter(X1 > 17, X1 < 50, X2 < -15, X2 > -50)
retina_plus_group <- bind_rows(retina_plus_group, tsne_50 %>% filter(X1 > 35))
dbscan_cluster <- dbscan(retina_plus_group %>%  select(X1, X2) ,  minPts = 3, eps = 2)

retina_plus_group$Cluster <- dbscan_cluster$cluster
col1 = 'X1'
col2 = 'X2'
# create label points for each cluster
cluster_centers <- retina_plus_group %>%
  left_join(.,core_tight_2019, by = c("sample_accession", "study_accession", "study_title", "study_abstract", "sample_attribute", "Tissue", "Sub_Tissue", "Origin", "Kept")) %>% group_by(Cluster) %>%
  summarise(C1=mean(!!sym(col1)),C2=mean(!!sym(col2)),Tissue=paste(unique(Sub_Tissue),collapse=','))

# samples closest to center in each cluster
center_samples <- retina_plus_group %>% left_join(.,core_tight_2019, by = c("sample_accession", "study_accession", "study_title", "study_abstract", "sample_attribute", "Tissue", "Sub_Tissue", "Origin", "Kept"))  %>%
  left_join(.,cluster_centers, by=c('Cluster')) %>%
  mutate(Distance = (!!sym(col1)-C1)+(!!sym(col2)-C2)) %>%
  group_by(Cluster) %>%
  dplyr::slice(which.min(Distance)) %>%
  dplyr::slice(1) %>%
  .[['sample_accession']]

# cluster stats
cluster_stats <- retina_plus_group %>% left_join(.,core_tight_2019, by = c("sample_accession", "study_accession", "study_title", "study_abstract", "sample_attribute", "Tissue", "Sub_Tissue", "Origin", "Kept"))  %>%
  mutate(Cluster = as.factor(Cluster)) %>%
  group_by(Cluster) %>%
  summarise(Cluster_Tissues = paste(unique(Sub_Tissue), collapse=',\n'), Cluster_Counts = paste(n(), ' samples', sep='')) %>% 
  mutate(Cluster_Tissues = paste0(Cluster, ': ', Cluster_Tissues))

# set up for ggplot
retina_plus_group_prep <- retina_plus_group %>%
  mutate(Origin=factor(Origin, levels=c('Adult Tissue', 'Fetal Tissue', 'Stem Cell', 'Cell Line', 'Organoid < 30 days', 'Organoid > 30 days'))) %>%
  mutate(Cluster = as.factor(Cluster)) %>%
  select(-Cluster_Tissues) %>%
  left_join(., cluster_stats, by=c('Cluster')) %>%
  mutate(Label = Cluster) %>%
  mutate(Label = ifelse(sample_accession %in% center_samples, Label, ""))




tsne_zoom <- retina_plus_group_prep  %>%
  mutate(Group = case_when(study_accession == 'SRP012682' ~ 'GTEx',
                           TRUE ~ Tissue),
         Label = case_when(grepl('^0', Label) ~ '', TRUE ~ Label),
         Group = factor(Group, levels = c('Cornea', 'ESC', 'Lens', 'Retina', 'RPE','GTEx'))) %>%
  mutate( zoom = TRUE )

# https://cran.r-project.org/web/packages/ggplot2/vignettes/extending-ggplot2.html
StatChull <- ggproto("StatChull", Stat,
                     compute_group = function(data, scales) {
                       data[chull(data$x, data$y), , drop = FALSE]
                     },
                     
                     required_aes = c("x", "y")
)
stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

zoom_plot <- tsne_50_prep %>% 
  rowwise() %>% 
  mutate(Age = suppressWarnings(as.integer(str_split(sample_attribute, '_')[[1]][2]))) %>% 
  mutate(Origin = as.character(Origin)) %>% 
  mutate(Origin = case_when(Age < 30 ~ 'Organoid < 30 days',
                            Age > 30 ~ 'Organoid > 30 days',
                            TRUE ~ Origin)) %>% 
  mutate(Origin=factor(Origin, levels=c('Adult Tissue', 'Fetal Tissue', 'Stem Cell', 'Cell Line', 'Organoid < 30 days', 'Organoid > 30 days')),
         Tissue = case_when(sample_accession == 'SRS1572426' ~ 'Retina',
                            TRUE ~ Tissue)) %>% 
  mutate(ZoomTissue = case_when((Tissue == 'Retina' & Sub_Tissue != 'Retina - Adult Tissue') | Tissue == 'ESC' ~ 'Early Retina',
                                TRUE ~ Sub_Tissue),
         Group = case_when(study_accession == 'SRP012682' ~ 'GTEx',
                           TRUE ~ Tissue),
         Label = case_when(grepl('^0', Label) ~ '', TRUE ~ Label)) %>% 
  mutate(Group = factor(Group, levels = c('Cornea', 'ESC', 'Lens', 'Retina', 'RPE','GTEx'))) %>% 
  ggplot(aes(x=X1,y=X2)) +
  scale_shape_manual(values=c(0:20,35:50)) +
  #geom_point(size=8, alpha=0.2, aes(colour=Group)) +
  geom_point(size=2, alpha=1, aes(shape=Origin, colour = Group)) +
  xlab('t-SNE 1') + ylab('t-SNE 2') +
  facet_zoom(xy = Group == 'Retina', zoom.size = 2, horizontal = T, zoom.data=zoom) + 
  geom_label_repel(data=tsne_zoom,size = 2, aes(label=Label), alpha=0.7, box.padding = unit(0.3, "lines"), force = 10) + 
  stat_chull(data=tsne_zoom %>% 
               mutate(`Sub-Tissue Cluster` = Cluster_Tissues),aes(fill=`Sub-Tissue Cluster`), alpha = 0.3) +
  scale_color_discrete_sequential(palette = 'viridis') + 
  scale_fill_viridis_d(option = 'plasma') +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(alpha = 1), ncol=2),
         shape = guide_legend(ncol =2 )) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        text = element_text(family = 'Linux Libertine O', size = 6))

heatmap <- cowplot::ggdraw() + cowplot::draw_image('figures_and_tables/heatmap_retina_time_series.svg')
svg("figures_and_tables/zoom_heatmap_retina.svg", width = 12, height = 12)
cowplot::plot_grid(zoom_plot, NULL, heatmap, ncol = 3, rel_widths = c(9,-2.3, 3), scale=c(0.5,0,1), margin)
dev.off()

