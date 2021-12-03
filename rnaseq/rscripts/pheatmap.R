BiocManager::install("pheatmap")
BiocManager::install("ComplexHeatmap")
library(pheatmap)
library(ComplexHeatmap)

setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/Circadian_Statistics")

circdata <- read.csv("Organized_by_time_genelist.csv", header = T, row.names = "Gene_ID")
#circdata <- read.csv("Organized_by_time_genelist.csv")

circdata <- as.matrix(circdata)

pheatmap(circdata, cluster_rows = TRUE, cluster_cols = FALSE)
pheatmap(circdata3, main = "Significantly Rhythmic Male Genes", scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE)

circdata2 <- read.csv("trial_graphdata.csv", header = T, row.names = "Gene_ID")
circdata2 <- as.matrix(circdata2)

pheatmap(circdata2, cluster_rows=FALSE, cluster_cols=FALSE)

circdata3 <- read.csv("trial_graphdata2.csv", header = T, row.names = "Gene_ID")
circdata3 <- as.matrix(circdata3)

pheatmap(circdata3, cluster_rows=TRUE, cluster_cols=FALSE)

mcounts <- read.csv("mcountsnu2.average.csv", header = T, row.names = "Gene_ID")
fcounts <- read.csv("fcountsnu2.average.csv", header = T, row.names = "Gene_ID")

mcounts <- as.matrix(mcounts)
fcounts <- as.matrix(fcounts)
mcounts2 <- read.csv("mcountsnu2.averagewmodules.csv", header = T, row.names = "Gene_ID")
mcounts.f <- as.matrix(mcounts2[,1:8])
module_membership <- as.matrix(mcounts2[,9:11])
subset_module <- as.matrix(mcounts2[,12])
module_membership2 <- module_membership[subset_module[,1]]
subset(module_membership[,2], subset_module)
organized_modules <- read.csv("Male_module_membership.csv")

mrow_order <- as.data.frame(organized_modules[,1])
module_names <- as.data.frame(organized_modules[,2])

row_order <- as.data.frame(mrow_order[,1])
rownames(row_order) = mrow_order[,1]
rownames(mcounts) = rownames(row_order)
match_test <- match(rownames(mcounts), rownames(row_order))
match_test <- as.data.frame(match_test)
modules = data.frame("Modules" = module_names)
names(modules) = "modules"
rownames(morganized) = rownames(mrow_order)

pheatmap(mcounts, main = "All Male Genes", scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_row = modules)

pheatmap(mcounts, main = "All Male Genes", scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE)
mh <- pheatmap(mcounts, main = "All Male Genes", scale = "row", cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE)

pheatmap(fcounts, main = "All Female Genes", scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE)
fh <- pheatmap(fcounts, main = "All Female Genes", scale = "row", cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE)
mh + fh
