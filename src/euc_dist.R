library(viridis)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(pool)
library(RSQLite)
gene_pool_2019 <- dbPool(drv = SQLite(), dbname = '~/git/eyeIntegration_app/www/2019/EiaD_human_expression_2019_04.sqlite')

core_tight_2019 <- gene_pool_2019 %>% tbl('metadata') %>% as_tibble()

gene <- all_genes
fetal_ret <- core_tight_2019 %>% filter(Sub_Tissue == 'Retina - Fetal Tissue') %>% pull(sample_accession)
query = paste0('select * from lsTPM_gene where sample_accession in ("',paste(fetal_ret, collapse='","'),'")')
p <- dbGetQuery(gene_pool_2019, query) %>% left_join(.,core_tight_2019) %>% 
  left_join(., gene_pool_2019 %>% tbl('gene_IDs') %>% as_tibble()) %>% 
  as_tibble()

x <- p %>% 
  group_by(ID, Age_Days) %>% 
  summarise(value = mean(value)) %>% 
  mutate(Days = as.integer(Age_Days)) %>% 
  select(-Age_Days) %>% 
  spread(ID, value) %>% t()
colnames(x) <- x[1,]
x <- x[-1,]

euc_dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
e<-apply(x, 1, function(y) euc_dist(x['OTX2',],y))
e %>% sort() %>% head()
