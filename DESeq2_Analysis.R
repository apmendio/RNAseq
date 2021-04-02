Author = "Aron Judd P. Mendiola"

# Raw text files can be found in the cluster directory path below
# "/share/lasallelab/Aron/RNA_seq_analysis/Analysis_March252021/GeneCounts/ReverseCounts"

# setwd = set working directory; equivalent to the Linux "cd".
# the R equivalent to the Linux pwd is getwd() = get working directory.
setwd("/Users/aron/Desktop/LaSalle Lab/Analysis/Deseq2/ReverseCounts")

# load package DESeq2 (all functions)
library(DESeq2)

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
dir(pattern="counts.txt")

# save the results to a variable
files <- dir(pattern="counts.txt")

counts <- c()
for( i in seq_along(files) ){
  x <- read.table(file=files[i], sep="\t", header=F, as.is=T)
  counts <- cbind(counts, x[,2])
}

# subsetting data for male and female samples

males <- subset(se_star_matrix2$S)

# set the row names
rownames(counts) <- x[,1]
# set the column names based on input file names, with pattern removed
colnames(counts) <- sub("counts.txt","",files)

# Option 1 that reads in a matrix (we will not do it here):

# first read in the matrix
count_matrix <- read.delim("counts", header=T, sep="\t", row.names=1)

# then create the DESeq object
# countData is the matrix containing the counts
# sampletable is the sample sheet / metadata we created
# design is how we wish to model the data: what we want to measure here is the difference between the treatment times
se_star_matrix <- DESeqDataSetFromMatrix(countData = counts,
                                         colData = sampletable,
                                         design = ~ Timepoint)

# Number of genes before filtering:
nrow(se_star_matrix)

# Filter
se_star_matrix <- se_star_matrix[rowSums(counts(se_star_matrix)) > 10, ]


# Number of genes left after low-count filtering:
nrow(se_star_matrix)

# Run Deseq2
se_star_matrix2 <- DESeq(se_star_matrix)

# compute normalized counts (log2 transformed); + 1 is a count added to avoid errors during the log2 transformation: log2(0) gives an infinite number, but log2(1) is 0.
# normalized = TRUE: divide the counts by the size factors calculated by the DESeq function
norm_counts <- log2(counts(se_star_matrix2, normalized = TRUE)+1)

# Load the tximport package that we use to import Salmon counts
library(tximport)

# Read in the two-column data.frame linking transcript id (column 1) to gene id (column 2)
tx2gene <- read.table("tx2gene.gencode.vm25.csv", 
                      sep="\t",
                      header=F)

gene_list <- gsub("\\..*","", tx2gene)

# add the gene symbols
norm_counts_symbols <- merge(unique(tx2gene[,2:3]), data.frame(ID=rownames(norm_counts), norm_counts), by=1, all=T, verbose=T)

norm_counts_symbols2 <- merge((tx2gene[,2:3]), data.frame(ID=rownames(norm_counts), norm_counts), by=1, all=T, verbose=T)

# write normalized counts to text file
write.table(norm_counts_symbols, "normalized_counts.txt", quote=F, col.names=T, row.names=T, sep="\t")

# Try with the vst transformation
vsd <- vst(se_star_matrix2)

# load libraries pheatmap to create the heatmap plot
library(pheatmap)

# calculate between-sample distance matrix
sampleDistMatrix <- as.matrix(dist(t(assay(vsd))))

# create figure in PNG format
png("sample_distance_heatmap_star.png")
pheatmap(sampleDistMatrix)
# close PNG file after writing figure in it
dev.off() 

png("PCA_star.png")
plotPCA(object = vsd,
        intgroup = "Timepoint")
dev.off()

# check results names: depends on what was modeled. Here it was the "Timepoint"
resultsNames(se_star_matrix2)

# extract results for t25 vs t0
# contrast: the column from the metadata that is used for the grouping of the samples (Time), then the baseline (t0) and the group compared to the baseline (t25) -> results will be as "t25 vs t0"
de <- results(object = se_star_matrix2, 
              name="Timepoint_ZT6_vs_ZT0")

# processing the same results as above but including the log2FoldChange shrinkage
# useful for visualization and gene ranking
de_shrink <- lfcShrink(dds = se_star_matrix2,
                       coef="Timepoint_ZT6_vs_ZT0",
                       type="apeglm")

# check first rows of both results
head(de)
head(de_shrink)

# write normalized counts to text file
write.table(de_shrink, "ZT6_vs_ZT0.txt", quote=T, col.names=T, row.names=T, sep="\t")

# column 4 is the log2FoldChange, column 7 is the adjusted p-value (padj)
# keep all columns
awk '($7 < 0.05 && $4 > 0.5) || ($7 < 0.05 && $4 < -0.5) {print}' ZT6_vs_ZT0.txt > ZT6_vs_ZT0_padj0.05_log2fc0.5.txt

# extract only gene IDs (column 1)
cut -f1 ZT6_vs_ZT0_padj0.05_log2fc0.5.txt > ZT6_vs_ZT0_padj0.05_log2fc0.5_IDs.txt

# extract only gene symbols (column 2)
cut -f2 ZT6_vs_ZT0_padj0.05_log2fc0.5.txt > ZT6_vs_ZT0_padj0.05_log2fc0.5_symbols.txt

