library(viridis)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(pool)
library(RSQLite)
gene_pool_2019 <- dbPool(drv = SQLite(), dbname = '~/git/eyeIntegration_app/www/2019/EiaD_human_expression_2019_04.sqlite')

core_tight_2019 <- gene_pool_2019 %>% tbl('metadata') %>% as_tibble()

# gene list ------
rgc <- c('GAP43', 'POU4F1', 'ISL1', 'POU4F2','ATOH7','DLX2','SHH','DLX2')
progenitor <- c('VSX2','SOX2','SOX9','ASCL1','SFRP2','HES1','LHX2','PRTG','LGR5','ZIC1','DLL3','GLI1','FGF19','LIN28B')
# cone_rod <- c('NEUROD1','CRX','RORB','GUCA1B','GUCA1A','GUCY2D','PRPH2','RP1','RBP3','TULP1','AIPL1','RCVRN','GUCY2F','SLC24A1')
cone <- c('RXRB','THRB','RORA','GNAT2','ARR3','GNGT2','PDE6C','CNGA3','PDE6H','GNB3','GUCA1C','OPN1MW','OPN1SW','OPN1LW','GRK7')
rod <- c('NR2E3','NRL','MEF2C','ESRRB','CNGB1','GNAT1','GNGT1','GRK1','PDE6G','PDE6A','CNGA1','RHO','SAG','GNB1','PDE6B')
pr <- c('OTX2','RCVRN','AIPL1','NRL','CRX','PDE6B','NR2E3','ROM1','GNGT2','GNAT1','PDE6H','CNGB1','OPN1SW','GUCA1A','GNAT2','CNGA1','RHO','OPN1MW', cone, rod) %>% unique()
all_markers <- c(rgc, pr, progenitor, cone, rod) %>% unique()


# fetal retina
fetal_ret <- core_tight_2019 %>% filter(Sub_Tissue == 'Retina - Fetal Tissue') %>% pull(sample_accession)
query = paste0('select * from lsTPM_gene where sample_accession in ("',paste(fetal_ret, collapse='","'),'")')
p <- dbGetQuery(gene_pool_2019, query) %>% left_join(.,core_tight_2019) %>% 
  left_join(., gene_pool_2019 %>% tbl('gene_IDs') %>% as_tibble()) %>% 
  as_tibble()

x <- p %>%
  select(ID, sample_accession, value) %>% 
  unique() %>% 
  spread(ID, value) %>% 
  select(-sample_accession) %>% 
  t()
sample_accession <- p %>%
  select(ID, sample_accession, value) %>% 
  unique() %>% 
  spread(ID, value) %>% 
  pull(sample_accession)
colnames(x) <- sample_accession
x <- x %>% as.matrix()

euc_dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

top <- ''
#bottom <- ''
for (i in all_markers){
  print(i)
  distances <- apply(x, 1, function(y) euc_dist(x[i,],y))
  distances <- distances %>% sort() %>% head(100)
  distances <- distances[!names(distances) %in% all_markers]
  head <- names(distances)[1]
  #tail <- distances %>% sort() %>% tail(1)
  top <- c(top, head)
  #tail <- c(bottom, tail)
}

link <- cbind(all_markers, top[2:length(top)])
colnames(link) <- c('Marker', 'TopRelated')
save(link, file = 'data/markers_related_retina.Rdata')










# ALL (for HPC)
library(parallelDist)
library(reshape)
query = paste0('select * from lsTPM_gene')
p <- dbGetQuery(gene_pool_2019, query) %>% left_join(.,core_tight_2019) %>% 
  left_join(., gene_pool_2019 %>% tbl('gene_IDs') %>% as_tibble()) %>% 
  as_tibble()
x <- p %>%
  select(ID, sample_accession, value) %>% 
  unique() %>% 
  spread(ID, value) %>% 
  select(-sample_accession) %>% 
  t()
sample_accession <- p %>%
  select(ID, sample_accession, value) %>% 
  unique() %>% 
  spread(ID, value) %>% 
  pull(sample_accession)
colnames(x) <- sample_accession
y <- parDist(x, method = 'euclidean', threads = 16) %>% as.matrix()
euc_dist_all_by_all <- melt(y)
euc_dist_all_by_all_top100 <- euc_dist_all_by_all %>% group_by(Var2) %>% arrange(value) %>% top_n(-100)
