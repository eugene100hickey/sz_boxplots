#just two libraries needed here
library(tidyverse)
require(ABAData)

#input data from ABAData
data("dataset_adult")

#need these two csv files
#they come from the supplementary tables from Ripke(2014)
#generates list of SZ genes
data <- read_csv("nature13595-s3.csv")
sz_genes <- data$'eQTL gene'
data <- read_csv("nature13595-s3-blood.csv")
sz_genes_blood <- data$'eQTL gene'
all_sz_genes <- c(sz_genes, sz_genes_blood)

#filter just for the genes from SZ
#recast ABAData in wide format
#remove incomplete and duplicate rows
dataset_adult <- dataset_adult %>% filter(hgnc_symbol %in% all_sz_genes)
wide_allen <- spread(dataset_adult, structure, signal)
wide_allen <- wide_allen[complete.cases(wide_allen),]
duplicate_index <- duplicated(wide_allen$entrezgene)
wide_allen <- wide_allen[!duplicate_index,]

#does a kmeans clustering with 4 groups
km <- kmeans(wide_allen_sz[,5:418], centers = 4)
wide_allen_sz <- wide_allen_sz %>% 
  mutate(cluster = as.factor(km$cluster))
#recasts dataframe into long format
long_allen_sz <- wide_allen_sz %>% 
  gather(area, expression, -c(hgnc_symbol:age_category, cluster))

#pretty boxplot picture
long_allen_sz %>% 
  ggplot(aes(group = cluster, y = expression, colour = cluster)) + 
  geom_boxplot() + 
  theme_minimal() + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())