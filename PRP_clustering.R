## clustering analysis of dataset resulting from mass spectrometry proteomic analysis of equine platelet-rich plasma (PRP)
# investigating associations between donor phenotype, protein profile of PRP and its functional implication 


# load packages
library(tidyverse)
library(ggplot2)
library(DEP)
library(SummarizedExperiment)

# load protein intensity vs sample data (MAxQuant LFQ_intensity)
prp <- read.csv('C:/Users/aggiej/Documents/Proteomics/PRP_VPET/PRP_LFQ_filtered.csv')

# rename the first column
colnames(prp)[1] <- 'ID'

## generate SummarisedExperiment object
# create column to be used as condition in exp design
prp$cond <- paste(prp$Age, prp$Sex, prp$Injury, sep = '_')

# extract data for assay object and format
prp_a <- as.data.frame(t(prp %>% select(!c('cond', 'ID', 'Age', 'Sex', 'Injury'))))
colnames(prp_a) <- prp$ID

# create required columns 'name' and 'ID'
prp_a$name <- rownames(prp_a)
prp_a <- prp_a %>% relocate(name, .before = 'PRP1')
prp_a <- prp_a %>% mutate(ID = row_number())
prp_a <- prp_a %>% relocate(ID, .after = 'name')

# create experimental design
prp_design <- prp %>% select(c('cond', 'ID'))

# 'label' column must correspond with column names in assay data frame
# 'condition' and 'replicate' columns are also required - values combination cannot be identical between any two rows
colnames(prp_design)[2] <- 'label'
colnames(prp_design)[1] <- 'condition'
prp_design <- prp_design %>% 
  tibble::as_tibble() %>% 
  group_by(condition) %>% 
  dplyr::mutate(replicate = row_number())

# create unique identifiers based on name and ID - prerequisite, even if names are already unique
unique_names <- make_unique(prp_a, 'name', 'ID', delim = ';')

# create SummarisedExperiment
prp_se <- make_se(unique_names, 3:30, prp_design)

# Scale and variance stabilize
prp_norm <- normalize_vsn(prp_se)

#
#
#

# unsupervised clustering
# non-negative matrix factorisation (NMF)
library(NMF)
library(stringr)
library(readr)

# extract metadata for later
col_names <- colnames(prp_norm)
age <- as.numeric(parse_number(col_names))
sex <- word(col_names, 2 , sep = "_")
injury <- word(col_names, 3 ,sep = "_")

# extract the abundance data from se object - remove NA containing rows
prp_norm_no_na <-na.omit(as.data.frame(assay(prp_norm)))

# run non-smooth NMF with 2-6 clusters - 100 runs each initially
estim.r <- nmf(prp_norm_no_na, 2:6, nrun = 100, seed = 123, method = "lee")

#plot cluster stability statistics to help choose the k
plot(estim.r)

# k = 4
res <- nmf(prp_norm_no_na, 4, nrun=200, seed=123,method = "lee")
summary(res)

# plot the consensus clustering along side the metadata
# open a new device
dev.new()

# plot clustering
consensusmap(res, annCol = data.frame("age" = age, "injury" = injury, "sex"=sex), tracks = c("basis:", "consensus:", "silhouette:"))
dev.off()

# plot the meta-genes and meta-samples
coefmap(res)
basismap(res)

#
#
#

# PCA

library(factoextra)
library(patchwork)

pca <-prcomp(t(prp_norm_no_na), scale=TRUE, center = TRUE)
pca_scores <- as.data.frame(pca$x)
pca_load <- as.data.frame(pca$rotation)

plot <- ggplot(pca_scores, aes(x=PC1, y=PC2, colour = prp$Sex)) + geom_point() + geom_text(aes(label = prp$ID)) + ylim(-7, 10)
plot
