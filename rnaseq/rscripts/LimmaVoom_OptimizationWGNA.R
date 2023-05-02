library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)

BiocManager::install("limma")
BiocManager::install("Glimma")
BiocManager::install("edgeR")
BiocManager::install("Mus.musculus")
BiocManager::install("")

library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
library()

setwd("/Users/aron/Desktop/LaSalle Lab/Analysis/Limma/ReverseCounts")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/primate_wgcna/chr10")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/oran")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/clamsrw/rnaseq/de")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/clamsrw/rnaseq/de2")
# read in the sample sheet
# header = TRUE: the first row is the "header", i.e. it contains the column names.
# sep = "\t": the columns/fields are separated with tabs.
sampletable <- read.table("sample_sheet2.txt", header=T, sep="\t")
baboon_counts <-read.csv("baboon_counts.csv")
clams_rw_counts <- read.csv("clams-rw_counts.csv")
sampletable <- read.table("clams-rw_traitsUpdated4.txt", header=T, sep="\t")
clams_rw_counts <- read.csv("clams-rw_counts.csv")

clams_rw_counts <- read.csv("clams-rw_counts.txt", header=T, sep="\t")
sampletable <- read.csv("clams-rw_traits.csv")

sampletable <- read.csv("clams_male.traits.csv")
clams_male_counts <- read.csv("clams_male.counts.csv")
# add the sample names as row names (it is needed for some of the DESeq functions)
rownames(sampletable) <- sampletable$SampleID

# display the first 6 rows
head(sampletable)

# check the number of rows and the number of columns
nrow(sampletable) # if this is not 6, please raise your hand !
ncol(sampletable) # if this is not 4, also raise your hand !

# run this if counts files are not combined
# only return file names with a given pattern
dir(pattern="__counts.txt")
  
# save the results to a variable
files <- dir(pattern="__counts.txt")

counts <- c()
for( i in seq_along(files) ){
  x <- read.table(file=files[i], sep="\t", header=F, as.is=T)
  counts <- cbind(counts, x[,2])
}

counts <- exp_femdata2
df <- exp_femdata[!duplicated(exp_femdata$Gene_ID),]

any(duplicated(counts$Gene_ID))
counts <- df
rownames(df) <- df$Gene_ID
df <- df[-c(1)]
df2 <- subsetData(df)
rownames(counts) <- counts$Gene_ID
counts <- counts[-c(1)]

# run here if counts files have been combined in shell

counts <- baboon_counts
rownames(counts) <- baboon_counts$GeneID
counts <- counts[-c(1)]

counts <- clams_rw_counts
rownames(counts) <- clams_rw_counts$EnsemblID
counts <- counts[-c(1)]

counts <- clams_male_counts
rownames(counts) <- clams_rw_counts$EnsemblID
counts <- counts[-c(1)]

counts <- clams_rw_counts
rownames(counts) <- clams_rw_counts$EnsemblID
counts <- counts[-c(1)]

# set the row names
rownames(counts) <- x[,1]
# set the column names based on input file names, with pattern removed (if generated skip to line 54)
colnames(counts) <- sub("counts.txt","",files)
dim(counts)
head(counts)

counts <- luhmes

#Create Differential Gene Expression List Object

counts <- counts[-c(1)]

#here!!
d0 <- DGEList(counts)

#Read in Annotation

anno <- read.delim("ensembl_mm_100.tsv",as.is=T)

dim(anno)
head(anno)
tail(anno)
any(duplicated(anno$Gene.stable.ID))

#Derive experiment metadata from the sample names
clams_rw_cf <- read.csv("clams-rw_cf.csv")
clams_rw_cf <- read.csv("clams-cf.csv")
counts <- clams_rw_cf
rownames(counts) <- clams_rw_cf$EnsemblID
counts <- counts[-c(1)]
sampletable <- read.csv("clams-cf_traits.csv")

sample_names <- colnames(counts)
metadata <- sampletable
colnames(metadata) <- c("SampleID", "Genotype", "GenotypeScores", "Entrainment", "Sex", "SexScore", "Timepoint")
metadata

sample_names <- colnames(counts)
metadata <- read.table("sample_sheet2.txt", header=T, sep="\t")
colnames(metadata) <- c("SampleName", "FileName", "Timepoint", "Sex", "Genotype")
metadata


sample_names <- colnames(counts)
metadata <- read.table("clams-rw_traitsUpdated4.txt", header=T, sep="\t")
colnames(metadata) <- c("SampleID", "Experiment", "Entrainment", "Sex", "Genotype")
metadata

sample_names <- colnames(counts)
metadata <- sampletable
metadata

#Derive experiment metadata from the sample names
sample_names <- colnames(counts)
metadata <- read.table("sample_info.txt", header=T, sep="\t")
colnames(metadata) <- c("SampleID", "Treatment", "Cell")
metadata

#CF
#Create new variable grouping combining group of timepoints
metadata$group <- interaction(metadata$Genotype, metadata$Entrainment)
table(metadata$group)
table(metadata$Genotype)
table(metadata$Entrainment)

#Create new variable grouping combining group of timepoints
metadata$group <- interaction(metadata$Genotype, metadata$Entrainment)
table(metadata$group)
table(metadata$Sex)
table(metadata$Timepoint)

#Create new variable grouping combining group of timepoints
metadata$group <- interaction(metadata$Treatment, metadata$Cell)
table(metadata$group)
table(metadata$Treatment)
table(metadata$Cell)

#Create new variable grouping combining group of timepoints
metadata$group <- interaction(metadata$Genotype, metadata$Entrainment)
table(metadata$group)
table(metadata$Entrainment)

#Normalization factor calculation which doesn't normalize data

d0 <- calcNormFactors(d0)
d0$samples
dim(d0)
#Filtering genes

cutoff <- 2
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(counts) # number of genes before cleanup
dim(d) # number of genes left

plotMDS(d, col = as.numeric(metadata$group), cex=1)
plotMDS(d, col = as.numeric(metadata$Genotype), cex=1)
plotMDS(d, col = as.numeric(metadata$group), cex=1)

plotMDS(d, col = as.numeric(metadata$group), cex=1)
plotMDS(d, col = as.numeric(metadata$Treatment), cex=1)
plotMDS(d, col = as.numeric(metadata$group), cex=1)

plotMDS(d, col = as.numeric(metadata$group), cex=1)
plotMDS(d, col = as.numeric(metadata$Experiment), cex=1)
plotMDS(d, labels = metadata$GenotypeScores, col = as.numeric(metadata$Genotype), cex=1)
col.geno <- metadata$Genotype
levels(col.geno) <- brewer.pal(nlevels(col.geno), "Set1")
plotMDS(d, labels = metadata$Genotype, col = col.geno, cex=1)
title(main = "Genotype Groups")
plotMDS(d, labels = metadata$Genotype, col = as.numeric(metadata$Entrainment), cex=1)
plotMDS(d, col = as.numeric(metadata$Sex), cex=1)
plotMDS(d, col = as.numeric(metadata$group), cex=1)
#Extracting "normalized" expression table
logcpm <- cpm(d, prior.count=2, log=TRUE)
write.table(logcpm,"clams-rw_normalized.txt",sep="\t",quote=F)
write.table(logcpm,"counts_normalizednu.txt",sep="\t",quote=F)
write.table(logcpm,"luhmes_counts.txt",sep="\t",quote=F)
logcpm <- data.frame(logcpm)
logcpm$genes <- row.names(logcpm)
logcpm$genes
cutoff <- 2
drop <- which(apply(cpm(counts), 1, max) < cutoff)
normcounts <- counts[-drop,]
dim(normcounts) # number of genes left
dim(counts)
logcpm <- cpm(normcounts, prior.count=2, log=TRUE)
papdata <- logcpm
write.table(logcpm,"baboon_normcounts.txt",sep="\t",quote=F)
baboon_normcounts = t(logcpm)

logcpm <- cpm(normcounts, prior.count=2, log=TRUE)
clamsrwdata <- logcpm
write.table(logcpm,"clamsrw.txt",sep="\t",quote=F)
baboon_normcounts = t(logcpm)
#Transpose counts for WGCNA and Circadian Analysis
input_mat = t(logcpm)
dim(input_mat)
organized_colnames <- read.table("organized_samples3.txt", header=F, sep="\t")  #formatting samples to keep table organized according to replicates
allcol_order <- as.data.frame(organized_colnames[,1])
row_order <- as.data.frame(allcol_order[,1])
input_mat <- input_mat[allcol_order[,1],]
write.csv(input_mat, file = "input_wgcna.csv", append = TRUE, quote = FALSE, sep = "\t")
#Voom transformation and calculation of variance weights
group <- metadata$group
treatment <- metadata$Genotype
treatment <- metadata$Entrainment
mm <- model.matrix(~0 + group)
head(mm)
mm2 <- model.matrix(~0 + group + treatment2)
head(mm2)
#Voom
d <- 
y <- voom(d, mm, plot = T)
y2 <- voom(d, mm2, plot = T)
#Fitting linear models in Limma
fit <- lmFit(y, mm)
efit <- eBayes(fit)
plotSA(efit, main = "Test")
head(coef(fit))


contr.matrix <- makeContrasts(
  WT12vDel12 = groupwt.wt.12 - grouphet.wt.12,
  WT12vTg12 = groupwt.wt.12 - groupwt.tg.12,
  WT12vDelTg12 = groupwt.wt.12 - grouphet.tg.12,
  levels = colnames(mm))
contr.matrix
fit <- lmFit(y, mm)
vfit <- contrasts.fit(fit, contrasts = contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main = "Test")
head(coef(fit))

summary(decideTests(fit))
#Specify which groups to compare using contrasts of timepoints
colnames(coef(fit))
contr <- makeContrasts(groupCM.0.0.M - groupCM.1.0.M, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp
tmp <- eBayes(tmp)
tmp

#Specify which groups to compare using contrasts of timepoints
contr <- makeContrasts(groupFemale.ZT0 - groupMale.ZT0, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp
tmp <- eBayes(tmp)
tmp

#Specify which groups to compare using contrasts of timepoints
contr <- makeContrasts(groupUndifferentiated.luhmes - groupDifferentiated.luhmes, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp
tmp <- eBayes(tmp)
tmp

summary(decideTests(fit))
tfit <- treat(tmp, lfc=1)
dt <- decideTests(tfit)
summary(dt)

#venn diagram
contrast.matrix <- makeContrasts(groupDifferentiated.luhmes, groupUndifferentiated.luhmes + groupDifferentiated.luhmes, levels = colnames(coef(fit)))
contrast.matrix
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
top.test <- topTable(fit2, coef = 1, sort.by = "P", n = 11834)
reults <- decideTests(fit2)
summary(reults)
head(reults)
head(top.table$Gene[de.common], n = 20)
de.common <- which(reults[,1]!=0 & reults[,2]!=0)
length(de.common)
vennDiagram(reults, names = c("Differentiated", "Undifferentiated"), main = "DE Genes Between Differentiated and Undifferentiated LUHMES", circle.col=c("blue", "red"), cex.main = 0.8)
write.fit(fit2, reults, file="results.txt")

#MultipleTestingAdjustment
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))
head(top.table, 50)
row.names(top.table)
#Annotation and adding in cpms 
top.table$Gene <- rownames(top.table)
table(top.table$Gene)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table,anno[match(top.table$Gene,anno$Gene.stable.ID),],logcpm[match(top.table$Gene,rownames(logcpm)),])

head(top.table)

write.table(top.table, file = "wtwtvhetwt.txt", row.names = F, sep = "\t", quote = F)

#Enhanced Volcano
EnhancedVolcano(top.table, lab = top.table$Gene, x = 'logFC', y = 'adj.P.Val')

#Heatmap
genes_of_interest <- read.csv("genes_of_interest.csv")
rownames(genes_of_interest) <- genes_of_interest$Gene
gene_list <- genes_of_interest$Gene
i <- which(logcpm$genes %in% gene_list)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(i, scale="row",
          labRow=i, labCol=metadata$group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

logcpm <- logcpm[-c(7)]

heatmap.2(as.matrix(logcpm[rownames(genes_of_interest),]), col=mycol, scale="row", trace="none", density.info="none", 
          margin=c(12,10), lhei=c(2,10), dendrogram="column")

top100 <- head(top.table, 100)
top
rownames(top100)
logcpm[rownames(top100),]

heatmap.2(as.matrix(logcpm[rownames(top100),]), col=mycol, scale="row", trace="none", density.info="none", 
          margin=c(12,10), lhei=c(2,10), dendrogram="column")

#GO Terms
top100$Gene %>%
  #dplyr::select(x) %>%
  #purrr::flatten() %>%
  enrichR::enrichr(c("GO_Biological_Process_2021",
                     "GO_Cellular_Component_2021",
                     "GO_Molecular_Function_2021",
                     "KEGG_2021_Human",
                     "Panther_2016",
                     "Reactome_2016")) %>% 
  #purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = "")) %T>%
  #purrr::map(~ dplyr::filter(., Adjusted.P.value < 0.05)) %>% 
  #purrr::map(~ dplyr::filter(., stringr::str_detect(Genes, ";"))) %>% 
  openxlsx::write.xlsx(file = "Top100DEG_enrichr.xlsx") %>%
  DMRichR::slimGO(tool = "enrichR",
                  annoDb = "org.Hs.eg.db",
                  plots = FALSE) %T>%
  openxlsx::write.xlsx(file = glue::glue("Top100DEG_rrvgo_enrichr.xlsx")) %>%
  DMRichR::GOplot() %>%
  ggplot2::ggsave("Top100DEG_enrichr_plot.pdf",
                  plot = .,
                  device = NULL,
                  height = 8.5,
                  width = 10) 
#Read differential expression data

#Above worked now test for multiple sets of analysis

#ZTO comparisons with all timepoints

contr <- makeContrasts(groupMale.ZT0 - groupMale.ZT3, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))
head(top.table, 20)
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table,anno[match(top.table$Gene,anno$Gene.stable.ID),],logcpm[match(top.table$Gene,rownames(logcpm)),])
head(top.table)
write.table(top.table, file = "ZT0_v_ZT3_Malenu.txt", row.names = F, sep = "\t", quote = F)

contr <- makeContrasts(groupMale.ZT0 - groupMale.ZT6, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))
head(top.table, 20)
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table,anno[match(top.table$Gene,anno$Gene.stable.ID),],logcpm[match(top.table$Gene,rownames(logcpm)),])
head(top.table)
write.table(top.table, file = "ZT0_v_ZT6_Malenu.txt", row.names = F, sep = "\t", quote = F)

contr <- makeContrasts(groupMale.ZT0 - groupMale.ZT9, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))
head(top.table, 20)
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table,anno[match(top.table$Gene,anno$Gene.stable.ID),],logcpm[match(top.table$Gene,rownames(logcpm)),])
head(top.table)
write.table(top.table, file = "ZT0_v_ZT9_Malenu.txt", row.names = F, sep = "\t", quote = F)

contr <- makeContrasts(groupMale.ZT0 - groupMale.ZT12, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))
head(top.table, 20)
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table,anno[match(top.table$Gene,anno$Gene.stable.ID),],logcpm[match(top.table$Gene,rownames(logcpm)),])
head(top.table)
write.table(top.table, file = "ZT0_v_ZT12_Malenu.txt", row.names = F, sep = "\t", quote = F)

contr <- makeContrasts(groupMale.ZT0 - groupMale.ZT15, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))
head(top.table, 20)
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table,anno[match(top.table$Gene,anno$Gene.stable.ID),],logcpm[match(top.table$Gene,rownames(logcpm)),])
head(top.table)
write.table(top.table, file = "ZT0_v_ZT15_Malenu.txt", row.names = F, sep = "\t", quote = F)

contr <- makeContrasts(groupMale.ZT0 - groupMale.ZT18, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))
head(top.table, 20)
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table,anno[match(top.table$Gene,anno$Gene.stable.ID),],logcpm[match(top.table$Gene,rownames(logcpm)),])
head(top.table)
write.table(top.table, file = "ZT0_v_ZT18_Malenu.txt", row.names = F, sep = "\t", quote = F)

contr <- makeContrasts(groupMale.ZT0 - groupMale.ZT21, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))
head(top.table, 20)
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table,anno[match(top.table$Gene,anno$Gene.stable.ID),],logcpm[match(top.table$Gene,rownames(logcpm)),])
head(top.table)
write.table(top.table, file = "ZT0_v_ZT21_Malenu.txt", row.names = F, sep = "\t", quote = F)

#ZT6 comparisons with all timepoints

contr <- makeContrasts(groupMale.ZT6 - groupMale.ZT12, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))
head(top.table, 20)
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table,anno[match(top.table$Gene,anno$Gene.stable.ID),],logcpm[match(top.table$Gene,rownames(logcpm)),])
head(top.table)
write.table(top.table, file = "ZT6_v_ZT12_Malenu.txt", row.names = F, sep = "\t", quote = F)

contr <- makeContrasts(groupMale.ZT6 - groupMale.ZT15, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))
head(top.table, 20)
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table,anno[match(top.table$Gene,anno$Gene.stable.ID),],logcpm[match(top.table$Gene,rownames(logcpm)),])
head(top.table)
write.table(top.table, file = "ZT6_v_ZT15_Malenu.txt", row.names = F, sep = "\t", quote = F)

contr <- makeContrasts(groupMale.ZT6 - groupMale.ZT21, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))
head(top.table, 20)
top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table,anno[match(top.table$Gene,anno$Gene.stable.ID),],logcpm[match(top.table$Gene,rownames(logcpm)),])
head(top.table)
write.table(top.table, file = "ZT6_v_ZT21_Malenu.txt", row.names = F, sep = "\t", quote = F)

#luhmes
library(ggplot2)
library(gplots)
top.table
signif.topgenes <- top.table$Gene[1:11834]
i <- which(top.table$Gene %in% signif.topgenes)
mycol <- colorpanel(1000, "blue", "white", "red")
heatmap.2(lcpm[i,], scale="row", labRow=top.table$Gene, labCol = group, col=mycol, trace ="none", density.info="none", margin=c(8,6), lhei=c(2,10), dendrogram="column")
