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

mcounts4 <- read.csv("mcountsnu2.average.csv", header = T, row.names = "Gene_ID")
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
rownames(modules) = rownames(mcounts4)
rownames(morganized) = rownames(mrow_order)

pheatmap(mcounts, main = "All Male Genes", scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_row = modules)

pheatmap(mcounts, main = "All Male Genes", scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE)
mh <- pheatmap(mcounts, main = "All Male Genes", scale = "row", cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE)

pheatmap(fcounts, main = "All Female Genes", scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE)
fh <- pheatmap(fcounts, main = "All Female Genes", scale = "row", cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE)
mh + fh

top.5.modules <- read.csv("top.5.modules.csv")
top.5.modules.mem <- read.csv("top.5.modules.mem.csv")
top.5.modules.mem <- as.data.frame(top.5.modules.mem[,2])
names(top.5.modules.mem) <- "Module Membership"
rownames(top.5.modules.mem) <- colnames(top.5.modules.df[,-c(1)])
top.5.modules.df <- as.data.frame(top.5.modules[,2:9])
rownames(top.5.modules.df) <- top.5.modules$Gene_ID
rownames(top.5.modules.mem) <- top.5.modules$Gene_ID

top.5.modules.mem <- as.data.frame(top.5.modules[,10:11])
pheatmap(top.5.modules.df, main ="Top 5 Male Modules", scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_row = top.5.modules.mem)

top.10.modules <- read.csv("modules.5to10.csv")
top.10.modules.mem <- as.data.frame(top.10.modules[,10:11])
names(top.10.modules.mem) <- "Module Membership"
top.10.modules.df <- as.data.frame(top.10.modules[,2:9])
rownames(top.10.modules.df) <- top.10.modules$Gene_ID
rownames(top.10.modules.mem) <- top.10.modules$Gene_ID

pheatmap(top.10.modules.df, main ="Top 10 Male Modules", scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_row = top.10.modules.mem)

top.15.modules <- read.csv("modules.11to15.csv")
top.15.modules.mem <- as.data.frame(top.15.modules[,10:11])
names(top.15.modules.mem) <- "Module Membership"
top.15.modules.df <- as.data.frame(top.15.modules[,2:9])
rownames(top.15.modules.df) <- top.15.modules$Gene_ID
rownames(top.15.modules.mem) <- top.15.modules$Gene_ID

top.15.modules.mem <- as.data.frame(top.15.modules[,10:11])
pheatmap(top.15.modules.df, main ="Top 15 Male Modules", scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_row = top.15.modules.mem)

top.20.modules <- read.csv("modules.16to20.csv")
top.20.modules.mem <- read.csv("top.20.modules.mem.csv")
top.20.modules.mem <- as.data.frame(top.20.modules.mem[,2])
names(top.20.modules.mem) <- "Module Membership"
rownames(top.20.modules.mem) <- colnames(top.20.modules.df[,-c(1)])
top.20.modules.df <- as.data.frame(top.20.modules[,2:9])
rownames(top.20.modules.df) <- top.20.modules$Gene_ID
rownames(top.20.modules.mem) <- top.20.modules$Gene_ID

top.20.modules.mem <- as.data.frame(top.20.modules[,10:11])
pheatmap(top.20.modules.df, main ="Top 20 Male Modules", scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_row = top.20.modules.mem)

top.25.modules <- read.csv("modules.21to25.csv")
top.25.modules.mem <- read.csv("top.25.modules.mem.csv")
top.25.modules.mem <- as.data.frame(top.25.modules.mem[,2])
names(top.25.modules.mem) <- "Module Membership"
rownames(top.25.modules.mem) <- colnames(top.25.modules.df[,-c(1)])
top.25.modules.df <- as.data.frame(top.25.modules[,2:9])
rownames(top.25.modules.df) <- top.25.modules$Gene_ID
rownames(top.25.modules.mem) <- top.25.modules$Gene_ID

top.25.modules.mem <- as.data.frame(top.25.modules[,10:11])
pheatmap(top.25.modules.df, main ="Top 25 Male Modules", scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, annotation_row = top.25.modules.mem)

