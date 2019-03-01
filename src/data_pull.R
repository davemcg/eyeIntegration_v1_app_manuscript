# pull data from pools and save as Rdata
library(pool)
library(RSQLite)
library(tidyverse)

gene_pool_2019 <- dbPool(drv = SQLite(), dbname = '~/git/eyeIntegration_app/inst/app/www/2019/EiaD_human_expression_2019_03.sqlite')
gene_pool_2017 <- dbPool(drv = SQLite(), dbname = '~/git/eyeIntegration_app/inst/app/www/2017/eyeIntegration_human_2017_01.sqlite')


core_tight_2017 <- gene_pool_2017 %>% tbl('metadata') %>% as_tibble()
core_tight_2019 <- gene_pool_2019 %>% tbl('metadata') %>% as_tibble()
load('data/human_tx_studies.Rdata')
tsne_50 <- gene_pool_2019 %>% tbl('tSNE_bulk_RNA') %>% filter(Perplexity == 45) %>% as_tibble()
gene_anno <- gene_pool_2019 %>% tbl('gene_IDs') %>% as_tibble()
tx_anno <- gene_pool_2019 %>% tbl('tx_IDs') %>% as_tibble()
limma_DE_tests <- gene_pool_2019 %>% tbl('limma_DE_tests') %>% as_tibble()
q01_limma_gene_sig <- gene_pool_2019 %>% tbl('limma_DE_gene') %>% filter(adj.P.Val < 0.01) %>% group_by(Comparison) %>% summarise(Count=n()) %>% as_tibble()
organoid_fetal_GO <- gene_pool_2019 %>% tbl('all_vs_all_GO') %>% filter(Set == 'Retina_Fetal.Tissue-Retina_3D.Organoid.Stem.Cell', ONTOLOGY == 'BP') %>% as_tibble() %>% mutate(qvalue = as.numeric(qvalue)) %>% arrange(qvalue)
save(core_tight_2017, 
     core_tight_2019, 
     tsne_50, gene_anno, 
     tx_anno, 
     limma_DE_tests, 
     q01_limma_gene_sig,
     organoid_fetal_GO,
     file = 'data/data_pull.Rdata')
