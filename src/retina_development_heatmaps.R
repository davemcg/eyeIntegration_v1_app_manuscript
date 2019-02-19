library(viridis)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(pool)
library(RSQLite)
gene_pool_2019 <- dbPool(drv = SQLite(), dbname = '~/git/eyeIntegration_app/www/2019/EiaD_human_expression_2019_04.sqlite')
rgc <- c('GAP43', 'POU4F1', 'ISL1', 'POU4F2','ATOH7','DLX2','SHH','DLX2')
pr <- c('OTX2','RCVRN','AIPL1','NRL','CRX','PDE6B','NR2E3','ROM1','GNGT2','GNAT1','PDE6H','CNGB1','OPN1SW','GUCA1A','GNAT2','CNGA1','RHO','OPN1MW')
progenitor <- c('VSX2','SOX2','SOX9','ASCL1','SFRP2','HES1','LHX2','PRTG','LGR5','ZIC1','DLL3','GLI1','FGF19','LIN28B')
cone_rod <- c('NEUROD1','CRX','RORB','GUCA1B','GUCA1A','GUCY2D','PRPH2','RP1','RBP3','TULP1','AIPL1','RCVRN','GUCY2F','SLC24A1')
cone <- c('RXRB','THRB','RORA','GNAT2','ARR3','GNGT2','PDE6C','CNGA3','PDE6H','GNB3','GUCA1C','OPN1MW','OPN1SW','OPN1LW','GRK7')
rod <- c('NR2E3','NRL','MEF2C','ESRRB','CNGB1','GNAT1','GNGT1','GRK1','PDE6G','PDE6A','CNGA1','RHO','SAG','GNB1','PDE6B')


all_markers <- c(rgc, pr, progenitor, cone_rod, cone, rod) %>% unique()

core_tight_2019 <- gene_pool_2019 %>% tbl('metadata') %>% as_tibble()

gene <- cone
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
               name = 'log2(TPM+1)',
               clustering_distance_rows = "pearson", 
               clustering_distance_columns = "euclidean")

x <- rbind(organoid_swaroop_GFP, ESC)
y <- x %>% select(-Type) %>% spread(ID, value) %>% t()
colnames(y) <- y['Days',]
colnames(y)[1] <- 'ESC'
y <- y[-1,]

two <- Heatmap(log2(y), cluster_columns = F, column_title = 'Kaewkhaw GFP+\n3D Retina',
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

three <- Heatmap(log2(y), cluster_columns = F, column_title = 'Kaewkhaw GFP-\n3D Retina',
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

################ 
# all in one
################
tissue <- bind_rows(fetal_tissue %>% mutate(Type = 'Tissue'), 
                    adult_tissue %>% mutate(Type = 'Tissue'),
                    ESC %>% mutate(Type = 'ESC'),
                    organoid_johnston %>% mutate(Type = 'Eldred'),
                    organoid_swaroop_GFP %>% mutate(Type = 'Kaewkhaw GFP+'),
                    organoid_swaroop_GFPneg %>% mutate(Type = 'Kaewkhaw GFP-'))
x <- tissue
y <- x %>% spread(ID, value) %>% select(-Days, -Type) %>% t()
type <- (x %>% spread(ID, value) %>% t())[2,]
days <- (x %>% spread(ID, value) %>% t())[1,]
colnames(y) <- days

ha_column = HeatmapAnnotation(df = data.frame(Type = type), 
                              col = list(Type = c("ESC" = magma(10)[1],
                                                  "Tissue" = magma(10)[3],
                                                  "Eldred" = magma(10)[5],
                                                  "Kaewkhaw GFP-" = magma(10)[7],
                                                  "Kaewkhaw GFP+" = magma(10)[9])))

Heatmap(log2(y), cluster_columns = T, 
        col = viridis(10),
        clustering_distance_rows = "pearson", 
        clustering_distance_columns = "euclidean", 
        top_annotation = ha_column,
        show_heatmap_legend = F)


# summarise by type -----
sum_type <- function(gene_set){
  query = paste0('select * from lsTPM_gene where ID in ("',paste(gene_set, collapse='","'),'")')
  p <- dbGetQuery(gene_pool_2019, query) %>% left_join(.,core_tight_2019) %>% 
    left_join(., gene_pool_2019 %>% tbl('gene_IDs') %>% as_tibble()) %>% 
    as_tibble()
  tissue <- p %>% 
    filter(Tissue == 'Retina' | Tissue == 'ESC') %>% 
    filter(Sub_Tissue != 'Retina - RGC Stem Cell') %>% 
    mutate(Age = case_when(Tissue == 'ESC' ~ 0,
                           Age_Days == '.' ~ 1000,
                           TRUE ~ as.numeric(Age_Days)),
           Type = case_when(study_accession == 'SRP159246' ~ 'Eldred',
                            grepl('GFP positive', sample_attribute) ~ 'Kaewkhaw GFP+',
                            grepl('GFP negative', sample_attribute) ~ 'Kaewkhaw GFP-',
                            Tissue == 'ESC' ~ 'ESC',
                            TRUE ~ 'Tissue')) %>% 
    group_by(ID, Age, Type) %>% 
    summarise(value = mean(value)) 
  
  
  x <- tissue# %>% filter(Type == type)
  y <- x %>% spread(ID, value) %>% t()
  age <- y['Age',]
  type <- y['Type',]
  z <-  x %>% spread(ID, value) %>% ungroup() %>% select(-Age, -Type) %>% t()
  colnames(z) <- age
  out <- list()
  out$data <- colSums(z)
  out$type <- type
  out
}

z <- rbind(
  sum_type(rgc)$data,
  sum_type(pr)$data,
  sum_type(progenitor)$data,
  sum_type(cone_rod)$data,
  sum_type(cone)$data,
  sum_type(rod)$data
) %>% as.matrix()
type <- sum_type(rgc)$type
row.names(z) <- c('RGC','PR','Progenitor','ConeRod','Cone','Rod')

ha_column = HeatmapAnnotation(df = data.frame(Type = type), 
                              col = list(Type = c("ESC" = magma(10)[1],
                                                  "Tissue" = magma(10)[3],
                                                  "Eldred" = magma(10)[5],
                                                  "Kaewkhaw GFP-" = magma(10)[7],
                                                  "Kaewkhaw GFP+" = magma(10)[9])))
Heatmap(log2(z %>% as.matrix()) , 
        cluster_columns = F, 
        col = viridis(10),
        top_annotation = ha_column,
        clustering_distance_rows = "pearson", 
        clustering_distance_columns = "euclidean", 
        show_heatmap_legend = F)


  # let's try merging ages into groups -----
ages <- colnames(z)
zz <- z %>% 
  t() %>% 
  as_tibble() 
zz$Age <- as.numeric(ages)
zz <- zz %>% mutate(Range = case_when(Age == 0 ~ 0,
                                Age < 40 ~ 35,
                                Age < 55 ~ 50,
                                Age < 60 ~ 60,
                                Age < 73 ~ 70,
                                Age == 80 ~ 80,
                                Age < 95 ~ 90,
                                Age < 106 ~ 100,
                                Age < 112 ~ 110,
                                Age < 120 ~ 120,
                                Age < 133 ~ 130,
                                Age < 145 ~ 140,
                                Age < 165 ~ 160,
                                Age < 185 ~ 180,
                                Age < 210 ~ 200,
                                Age < 230 ~ 220,
                                Age < 260 ~ 250,
                                TRUE ~ 1000))
zzz <- zz %>% group_by(Range) %>% summarise(RGC = mean(RGC), PR = mean(PR), Progenitor = mean(Progenitor), ConeRod = mean(ConeRod), Cone = mean(Cone), Rod = mean(Rod))
row.names(zzz) <- zzz$Range

Heatmap(log2(zzz %>% select(-Range) %>% as.matrix() %>% t()) , 
        cluster_columns = T, 
        col = viridis(10),
        #top_annotation = ha_column,
        clustering_distance_rows = "pearson", 
        clustering_distance_columns = "euclidean", 
        show_heatmap_legend = F)
######################





###################

# build model to predict fetal tissue age
# apply to organoid 
# doesn't really work

###################
tissue <- bind_rows(ESC, fetal_tissue, adult_tissue)

query = 'select * from lsTPM_gene'
p <- dbGetQuery(gene_pool_2019, query) %>% left_join(.,core_tight_2019) %>% 
  left_join(., gene_pool_2019 %>% tbl('gene_IDs') %>% as_tibble()) %>% 
  as_tibble()
tissue_allGene <- p %>% 
  filter(Sub_Tissue == 'Retina - Fetal Tissue') %>% 
  group_by(ID, Age_Days) %>% 
  summarise(value = mean(value)) %>% 
  mutate(Days = as.integer(Age_Days)) %>% 
  select(-Age_Days)
x <- tissue_allGene
y <- x %>% spread(ID, value) %>% t()
colnames(y) <- y['Days',]
y <- y[-1,]
# apply(y, 1, var) %>% enframe() %>% arrange(-value)
# top100 <- apply(y, 1, var) %>% enframe() %>% arrange(-value) %>% head(100) %>% pull(name)

data = y[all_genes,] %>% t() %>% data.frame()
#data
data$age <- row.names(data) %>% as.numeric()

# caret
library(caret)
library(xgboost)
set.seed(96)
lm_fit <- train(age ~ .,
                data = data, 
                method = "lm")

glmboost_fit <- train(age ~ .,
                      data = data, 
                      method = "glmboost")

rf_fit <- train(age ~ .,
                      data = data, 
                      method = "rf")
bst <- xgboost(data = data %>% select(-age) %>% as.matrix(), 
               label = data$age)
               , 
               max.depth = 3, 
               eta = 0.2, 
               gamma = 5,
               nrounds = 20,
               nthread = 2)


# predict on organoid

tissue <- p %>% 
  filter(study_accession == 'SRP159246') %>% 
  group_by(ID, Age_Days) %>% 
  summarise(value = mean(value)) %>% 
  mutate(Days = as.integer(Age_Days)) %>% 
  select(-Age_Days)
x <- tissue
y <- x %>% spread(ID, value) %>% t()
colnames(y) <- y['Days',]
y <- y[-1,]

predict(lm_fit, y %>% t())
predict(rf_fit, y %>% t())
predict(glmboost_fit, y %>% t())
predict(bst, y[colnames(data %>% select(-age)),] %>% t())
