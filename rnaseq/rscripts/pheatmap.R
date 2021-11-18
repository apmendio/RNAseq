BiocManager::install("pheatmap")
BiocManager::install("ComplexHeatmap")
library(pheatmap)
library(ComplexHeatmap)
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

pheatmap(mcounts, main = "All Male Genes", scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE)
mh <- pheatmap(mcounts, main = "All Male Genes", scale = "row", cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE)

pheatmap(fcounts, main = "All Female Genes", scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE)
fh <- pheatmap(fcounts, main = "All Female Genes", scale = "row", cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE)
mh + fh
