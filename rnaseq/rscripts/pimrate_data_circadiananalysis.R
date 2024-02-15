# Circadian Rhythm Analysis
#Load the package
library(BiocManager)
#Require rain package from biomanager
BiocManager::install("rain")
library(rain)
BiocManager::install("tidyverse")
library(tidyverse)
BiocManager::install("magrittr")
library(magrittr)
library(tidyverse)
library(oligo)
library(sva)
library(pd.hugene.2.0.st)
library(matrixStats)
library(biomaRt)
library(WGCNA)
library(VennDiagram)
BiocManager::install("tidyverse")
BiocManager::install("oligo")
BiocManager::install("sva")
BiocManager::install("pd.hugene.2.0.st")
BiocManager::install("matrixStats")
BiocManager::install("biomaRt")
BiocManager::install("WGCNA")
BiocManager::install("VennDiagram")
BiocManager::install("magrittr") # package installations are only needed the first time you use it
BiocManager::install("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
BiocManager::install("ggplot2")
library(ggplot2)
BiocManager::install("biomaRt")
library(biomaRt)

getwd()
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/primate_wgcna/test")
baboon_dataset <- read.csv("/Users/aron/Desktop/LaSalle_Lab/Analysis/primate_wgcna/test/baboon_prc_expression_FPKM.csv")
#change dataframe to matrix
test<-data.frame(baboon_dataset)
head(baboon_dataset)
colnames(baboon_dataset)
rownames(baboon_dataset)
rownames(test)<-baboon_dataset$EnsemblID
test2 <- test[-c(1,2)]
test2 <- t(test2)
#Use rain
results <- rain(test2, deltat=2, period=24, measure.sequence = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), peak.border=c(0.3, 0.7), verbose=FALSE)
#look at RAIN test result
results
BiocManager::install("biomaRt")
library(biomaRt)
listEnsembl()
# ensembl <- useEnsembl(biomart = "genes", dataset = "panubis_gene_ensembl")
ensembl <- useMart(biomart = "ensembl", dataset = "panubis_gene_ensembl")
rownames(results)
test3 <- results
test3$ensembl <- rownames(results)
genes_baboon <- lapply(test3, function(x){
  getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", values = x, mart = ensembl,
        verbose = TRUE) %>% unlist %>% as.character %>% unique %>% sort %>% paste(collapse = ", ")}) %>% unlist
#Print RAIN test result

write.csv(results, 'rhythmic_baboon.csv')

#

test3 <- read.csv("/Users/aron/Desktop/LaSalle_Lab/Analysis/primate_wgcna/test/rhythmic_baboon.csv")

ensembl <- useMart(biomart = "ensembl", dataset = "panubis_gene_ensembl")
rownames(results)
test3 <- results
test3$ensembl <- rownames(results)
genes_baboon <- lapply(test3, function(x){
  getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", values = x, mart = ensembl,
        verbose = TRUE) %>% unlist %>% as.character %>% unique %>% sort %>% paste(collapse = ", ")}) %>% unlist

