
library(tidyverse)

library(pool)
library(RSQLite)
gene_pool_2019 <- dbPool(drv = SQLite(), dbname = '/Volumes/Arges/eyeIntegration_app/www/2019/EiaD_human_expression_2019_03.sqlite')


###################

# build model to predict fetal tissue age
# apply to organoid 
# works OK on fetal tissue, but doesn't 
# work well with organoid

###################
tissue <- bind_rows(ESC, fetal_tissue, adult_tissue)

query = 'select * from lsTPM_gene'
p <- dbGetQuery(gene_pool_2019, query) %>% left_join(.,core_tight_2019) %>% 
  left_join(., gene_pool_2019 %>% tbl('gene_IDs') %>% as_tibble()) %>% 
  as_tibble()
fetalTissue_allGene <- p %>% 
  filter(gene_type == 'protein_coding') %>% 
  filter(Sub_Tissue == 'Retina - Fetal Tissue' | study_accession == 'SRP159246') %>% 
  select(ID, sample_accession, Age_Days, value) %>% 
  mutate(age_accession = paste0(Age_Days, sep='_', sample_accession)) %>% 
  select(-sample_accession, -Age_Days) %>% 
  unique()

# set up train test
age_access <- fetalTissue_allGene$age_accession %>% unique() %>% sort()
train_samp <- age_access[seq(1, length(age_access), 2)]
train <- fetalTissue_allGene %>% filter(age_accession %in% train_samp)
test <- fetalTissue_allGene %>% filter(!age_accession %in% train_samp)

# mean summarise multiple overlappping time points
train$age <- as.numeric(gsub('_.*','', train$age_accession))
train <- train %>% group_by(ID, age) %>% summarise(value = mean(value))

trainer <- train %>% spread(ID, value) 
tester <- test %>% spread(ID, value)

# caret
library(caret)
library(xgboost)
set.seed(96)

# rf_fit <- train(x = trainer %>% select(-age) %>% as.matrix(), 
#                 y = trainer$age,
#                 method = "rf")
xgb <- xgboost(data = trainer %>% select(-age) %>% as.matrix(), 
               label = trainer$age,
               nrounds = 50)

# predict on test
## rf
predictions <- predict(rf, tester %>% select(-age_accession) %>% as.matrix())
actual <- gsub('_.*','', tester$age_accession)
cbind(predictions, actual)

## xgboost
predictions <- predict(xgb, tester %>% select(-age_accession) %>% as.matrix())
actual <- gsub('_.*','', tester$age_accession)
cbind(predictions, actual) %>% as_tibble() %>% mutate_all(as.numeric) %>% ggplot(aes(x=actual, y=predictions)) + geom_point() + coord_cartesian(xlim=c(0,250), ylim=c(0,250)) + geom_abline(intercept = 0, slope = 1)+ geom_smooth(method = 'lm')
# predict on swaroop organoid
organoidTissue_allGene <- p %>% 
  filter(gene_type == 'protein_coding') %>% 
  filter(Sub_Tissue == 'Retina - 3D Organoid Stem Cell', study_accession != 'SRP159246') %>% 
  select(ID, sample_accession, Age_Days, value) %>% 
  mutate(age_accession = paste0(Age_Days, sep='_', sample_accession)) %>% 
  select(-sample_accession, -Age_Days) %>% 
  unique() %>% 
  spread(ID, value)

predictions <- predict(xgb, organoidTissue_allGene %>% select(-age_accession) %>% as.matrix())
actual <- gsub('_.*','', organoidTissue_allGene$age_accession)
cbind(predictions, actual)
cbind(predictions, actual) %>% as_tibble() %>% mutate_all(as.numeric) %>% ggplot(aes(x=actual, y=predictions)) + geom_point() + coord_cartesian(xlim=c(0,250), ylim=c(0,250)) + geom_abline(intercept = 0, slope = 1) + geom_smooth(method = 'lm')