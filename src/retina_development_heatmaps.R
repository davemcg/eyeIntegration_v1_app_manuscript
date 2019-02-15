library(viridis)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(pool)
library(RSQLite)
gene_pool_2019 <- dbPool(drv = SQLite(), dbname = '/Volumes/McGaughey_S/eyeIntegration_app/www/2019/EiaD_human_expression_2019_03.sqlite')
rgc <- c('GAP43', 'POU4F1', 'ISL1', 'POU4F2','ATOH7','DLX2','SHH','SLX2')
pr <- c('OTX2','RCVRN','AIPL1','NRL','CRX','PDE6B','NR2E3','ROM1','GNGT2','GNAT1','PDE6H','CNGB1','OPN1SW','GUCA1A','GNAT2','CNGA1','RHO','OPN1MW')
progenitor <- c('VSX2','SOX2','SOX9','ASCL1','SFRP2','HES1','LHX2','PRTG','LGR5','ZIC1','DLL3','GLI1','FGF19','LIN28B')

gene <- pr
query = paste0('select * from lsTPM_gene where ID in ("',paste(gene, collapse='","'),'")')
p <- dbGetQuery(gene_pool_2019, query) %>% left_join(.,core_tight_2019) %>% 
  left_join(., gene_pool_2019 %>% tbl('gene_IDs') %>% as_tibble()) %>% 
  as_tibble()

ESC <- p %>% 
  filter(Tissue == 'ESC') %>% 
  mutate(Days = 0, Type = 'ESC') %>% 
  group_by(ID, Days) %>% 
  summarise(value = mean(value)) %>% 
  mutate(Days = as.integer(Days))
organoid_swaroop_GFP <- p %>% 
  filter(Sub_Tissue == 'Retina - 3D Organoid Stem Cell', !grepl('GFP negative', sample_attribute), study_accession != 'SRP159246') %>% 
  group_by(ID, Age_Days) %>% 
  summarise(value = mean(value)) %>% 
  mutate(Days = as.integer(Age_Days), Type = 'GFP+ 3D Organoid') %>% 
  select(-Age_Days)
organoid_swaroop_GFPneg <-  p %>% 
  filter(Sub_Tissue == 'Retina - 3D Organoid Stem Cell', grepl('GFP negative', sample_attribute), study_accession != 'SRP159246') %>% 
  group_by(ID, Age_Days) %>% 
  summarise(value = mean(value)) %>% 
  mutate(Days = as.integer(Age_Days), Type = 'Kaewkhaw GFP- 3D Retina')%>% 
  select(-Age_Days)
organoid_johnston <-  p %>% 
  filter(study_accession == 'SRP159246') %>% 
  group_by(ID, Age_Days) %>% 
  summarise(value = mean(value)) %>% 
  mutate(Days = as.integer(Age_Days), Type = 'Kaewkhaw GFP+ 3D Retina') %>% 
  select(-Age_Days)
fetal_tissue <- p %>% 
  filter(Sub_Tissue == 'Retina - Fetal Tissue') %>% 
  group_by(ID, Age_Days) %>% 
  summarise(value = mean(value)) %>% 
  mutate(Days = as.integer(Age_Days), Type = 'Fetal Tissue') %>% 
  select(-Age_Days)
adult_tissue <- p %>% 
  filter(Sub_Tissue == 'Retina - Adult Tissue') %>% 
  group_by(ID) %>% 
  summarise(value = mean(value), Type = 'Adult Tissue') %>% 
  mutate(Days = 1000) 

tissue <- bind_rows(fetal_tissue, adult_tissue)
x <- tissue
y <- x %>% select(-Type) %>% spread(ID, value) %>% t()
colnames(y) <- y['Days',]
colnames(y)[ncol(y)] <- 'Adult'
y <- y[-1,]

one <- Heatmap(log2(y), cluster_columns = F,   column_title = 'Retina Tissue',
               col = viridis(10),
               show_row_names = FALSE,
               clustering_distance_rows = "pearson", 
               clustering_distance_columns = "euclidean")

x <- rbind(organoid_swaroop_GFP, ESC)
y <- x %>% select(-Type) %>% spread(ID, value) %>% t()
colnames(y) <- y['Days',]
colnames(y)[1] <- 'ESC'
y <- y[-1,]

two <- Heatmap(log2(y), cluster_columns = F, column_title = 'Kaewkhaw 3D\nOrganoid\nGFP Neg',
               col = viridis(10),
               clustering_distance_rows = "pearson", 
               clustering_distance_columns = "euclidean", 
               show_row_names = FALSE,
               show_heatmap_legend = F)

x <- rbind(organoid_swaroop_GFPneg, ESC)
y <- x %>% select(-Type) %>% spread(ID, value) %>% t()
colnames(y) <- y['Days',]
colnames(y)[1] <- 'ESC'
y <- y[-1,]

three <- Heatmap(log2(y), cluster_columns = F, column_title = 'Kaewkhaw 3D\nOrganoid\nGFP Neg',
                 col = viridis(10),
                 clustering_distance_rows = "pearson", 
                 clustering_distance_columns = "euclidean", 
                 show_row_names = FALSE,
                 show_heatmap_legend = F)

x <- rbind(organoid_johnston, ESC)
y <- x %>% select(-Type) %>% spread(ID, value) %>% t()
colnames(y) <- y['Days',]
colnames(y)[1] <- 'ESC'
y <- y[-1,]

four <- Heatmap(log2(y), cluster_columns = F, column_title = 'Eldred 3D Organoid', 
                col = viridis(10),
                clustering_distance_rows = "pearson", 
                clustering_distance_columns = "euclidean", 
                show_heatmap_legend = F)

one + two + three + four

