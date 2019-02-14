library(viridis)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
rgc <- c('GAP43', 'POU4F1', 'ISL1', 'POU4F2','ATOH7','DLX2','SHH','SLX2')
pr <- c('OTX2','RCVRN','AIPL1','NRL','CRX','PDE6B','NR2E3','ROM1','GNGT2','GNAT1','PDE6H','CNGB1','OPN1SW','GUCA1A','GNAT2','CNGA1','RHO','OPN1MW')
progenitor <- c('VSX2','SOX2','SOX9','ASCL1','SFRP2','HES1','LHX2','PRTG','LGR5','ZIC1','DLL3','GLI1','FGF19','LIN28B')

gene <- progenitor
query = paste0('select * from lsTPM_gene where ID in ("',paste(gene, collapse='","'),'")')
p <- dbGetQuery(gene_pool_2019, query) %>% left_join(.,core_tight_2019) %>% 
  left_join(., gene_pool_2019 %>% tbl('gene_IDs') %>% as.tibble()) %>% 
  as.tibble()

ESC <- p %>% 
  filter(Tissue == 'ESC') %>% 
  mutate(Days = 0) %>% group_by(ID, Days) %>% summarise(value = mean(value)) %>% mutate(Days = as.integer(Days))
organoid_swaroop_GFP <- p %>% 
  filter(sample_accession %in% (meta_ret %>% filter(Sub_Tissue == 'Retina - 3D Organoid Stem Cell') %>% pull(sample_accession))) %>%
  filter(!grepl('GFP negative', sample_attribute), study_accession != 'SRP159246') %>% 
  left_join(meta_ret, by = 'sample_accession') %>% group_by(ID, Days) %>% summarise(value = mean(value)) %>% mutate(Days = as.integer(Days))
organoid_swaroop_GFPneg <- p %>% 
  filter(sample_accession %in% (meta_ret %>% filter(Sub_Tissue == 'Retina - 3D Organoid Stem Cell') %>% pull(sample_accession))) %>%
  filter(grepl('GFP negative', sample_attribute), study_accession != 'SRP159246') %>% 
  left_join(meta_ret, by = 'sample_accession') %>% group_by(ID, Days) %>% summarise(value = mean(value)) %>% mutate(Days = as.integer(Days))
organoid_johnston <- p %>% 
  filter(sample_accession %in% (meta_ret %>% filter(Sub_Tissue == 'Retina - 3D Organoid Stem Cell') %>% pull(sample_accession))) %>%
  filter(!grepl('GFP negative', sample_attribute), study_accession == 'SRP159246') %>% 
  left_join(meta_ret, by = 'sample_accession') %>% group_by(ID, Days) %>% summarise(value = mean(value)) %>% mutate(Days = as.integer(Days))
  
tissue <- p %>% filter(sample_accession %in% (meta_ret %>% filter(Sub_Tissue == 'Retina - Fetal Tissue') %>% pull(sample_accession))) %>% 
  filter(!grepl('GFP negative', sample_attribute)) %>% 
  left_join(meta_ret, by = 'sample_accession') %>% group_by(ID, Days) %>% summarise(value = mean(value)) %>% mutate(Days = as.integer(Days))

x <- tissue
y <- x %>% spread(ID, value) %>% t()
colnames(y) <- y['Days',]
y <- y[-1,]

one <- Heatmap(log2(y), cluster_columns = F,   column_title = 'Retina Fetal Tissue',
               col = colorRamp2(c(-0, 15), viridis(2)))

x <- rbind(organoid_swaroop_GFP, ESC)
y <- x %>% spread(ID, value) %>% t()
colnames(y) <- y['Days',]
y <- y[-1,]

two <- Heatmap(log2(y), cluster_columns = F, column_title = 'Swaroop 3D Organoid\nGFP Neg',
               col = colorRamp2(c(-0, 15), viridis(2)),
               show_heatmap_legend = F)

x <- rbind(organoid_swaroop_GFPneg, ESC)
y <- x %>% spread(ID, value) %>% t()
colnames(y) <- y['Days',]
y <- y[-1,]

three <- Heatmap(log2(y), cluster_columns = F, column_title = 'Swaroop 3D Organoid\nGFP Neg',
               col = colorRamp2(c(-0, 15), viridis(2)),
               show_heatmap_legend = F)

x <- rbind(organoid_johnston, ESC)
y <- x %>% spread(ID, value) %>% t()
colnames(y) <- y['Days',]
y <- y[-1,]

four <- Heatmap(log2(y), cluster_columns = F, column_title = 'Johnston 3D Organoid', 
                 col = colorRamp2(c(0, 15), viridis(2)),
                show_heatmap_legend = F)

one + two + three + four
