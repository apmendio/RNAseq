library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)

BiocManager::install("limma")
BiocManager::install("Glimma")
BiocManager::install("edgeR")
BiocManager::install("Mus.musculus")

setwd("/Users/aron/Desktop/LaSalle Lab/Analysis/Limma/ReverseCounts")

# read in the sample sheet
# header = TRUE: the first row is the "header", i.e. it contains the column names.
# sep = "\t": the columns/fields are separated with tabs.
sampletable <- read.table("sample_sheet.txt", header=T, sep="\t")

# add the sample names as row names (it is needed for some of the DESeq functions)
rownames(sampletable) <- sampletable$SampleName

# display the first 6 rows
head(sampletable)

# check the number of rows and the number of columns
nrow(sampletable) # if this is not 6, please raise your hand !
ncol(sampletable) # if this is not 4, also raise your hand !

# only return file names with a given pattern
dir(pattern="__counts.txt")

# save the results to a variable
files <- dir(pattern="__counts.txt")

counts <- c()
for( i in seq_along(files) ){
  x <- read.table(file=files[i], sep="\t", header=F, as.is=T)
  counts <- cbind(counts, x[,2])
}

# set the row names
rownames(counts) <- x[,1]
# set the column names based on input file names, with pattern removed (if generated skip to line 54)
colnames(counts) <- sub("counts.txt","",files)
dim(counts)
head(counts)

#Create Differential Gene Expression List Object

d0 <- DGEList(counts)

#Read in Annotation

anno <- read.delim("ensembl_mm_100.tsv",as.is=T)

dim(anno)
head(anno)
tail(anno)
any(duplicated(anno$Gene.stable.ID))

#Derive experiment metadata from the sample names
sample_names <- colnames(counts)
metadata <- read.table("sample_sheet.txt", header=T, sep="\t")
colnames(metadata) <- c("SampleName", "FileName", "Timepoint", "Sex", "Genotype")
metadata

#Create new variable grouping combining group of timepoints
metadata$group <- interaction(metadata$Sex, metadata$Timepoint)
table(metadata$group)
table(metadata$Sex)
table(metadata$Timepoint)

#Normalization factor calculation which doesn't normalize data

d0 <- calcNormFactors(d0)
d0$samples

#Filtering genes

cutoff <- 2
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d) # number of genes left

plotMDS(d, col = as.numeric(metadata$group), cex=1)

#Extracting "normalized" expression table
logcpm <- cpm(d, prior.count=2, log=TRUE)
write.table(logcpm,"counts_normalized.txt",sep="\t",quote=F)

#Transpose counts for WGCNA and Circadian Analysis
input_mat = t(logcpm)
dim(input_mat)
organized_colnames <- read.table("organized_samples2.txt", header=F, sep="\t")  #formatting sampkes to keep table organized according to replciates
allcol_order <- as.data.frame(organized_colnames[,1])
row_order <- as.data.frame(allcol_order[,1])
input_mat <- input_mat[allcol_order[,1],]
#Voom transformation and calculation of variance weights
group <- metadata$group
sex <- metadata$Sex
mm <- model.matrix(~0 + group + sex)
head(mm)

#Voom
y <- voom(d, mm, plot = T)

#Fitting linear models in Limma
fit <- lmFit(y, mm)
head(coef(fit))

#Specify which groups to compare using contrasts of timepoints
contr <- makeContrasts(groupFemale.ZT0 - groupMale.ZT0, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)
tmp
tmp <- eBayes(tmp)
tmp

#MultipleTestingAdjustment
top.table <- topTable(tmp, adjust.method = "BH", sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))
head(top.table, 50)

#Annotation and adding in cpms 
top.table$Gene <- rownames(top.table)
table(top.table$Gene)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
top.table <- data.frame(top.table,anno[match(top.table$Gene,anno$Gene.stable.ID),],logcpm[match(top.table$Gene,rownames(logcpm)),])

head(top.table)

write.table(top.table, file = "ZT0_Female_v_Male.txt", row.names = F, sep = "\t", quote = F)
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
write.table(top.table, file = "ZT0_v_ZT3_Male.txt", row.names = F, sep = "\t", quote = F)

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
write.table(top.table, file = "ZT0_v_ZT6_Male.txt", row.names = F, sep = "\t", quote = F)

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
write.table(top.table, file = "ZT0_v_ZT9_Male.txt", row.names = F, sep = "\t", quote = F)

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
write.table(top.table, file = "ZT0_v_ZT12_Male.txt", row.names = F, sep = "\t", quote = F)

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
write.table(top.table, file = "ZT0_v_ZT15_Male.txt", row.names = F, sep = "\t", quote = F)

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
write.table(top.table, file = "ZT0_v_ZT18_Male.txt", row.names = F, sep = "\t", quote = F)

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
write.table(top.table, file = "ZT0_v_ZT21_Male.txt", row.names = F, sep = "\t", quote = F)

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
write.table(top.table, file = "ZT6_v_ZT12_Male.txt", row.names = F, sep = "\t", quote = F)

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
write.table(top.table, file = "ZT6_v_ZT15_Male.txt", row.names = F, sep = "\t", quote = F)

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
write.table(top.table, file = "ZT6_v_ZT21_Male.txt", row.names = F, sep = "\t", quote = F)
