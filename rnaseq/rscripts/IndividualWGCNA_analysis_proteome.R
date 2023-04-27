# Packages ####
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/individualWGCNA")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/primate_wgcna")
sapply(c("tidyverse", "oligo", "sva","pd.hugene.2.0.st", "matrixStats", "biomaRt", "WGCNA", "VennDiagram"), require, character.only = TRUE)

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
BiocManager::install("ggplot2")

library(tidyverse)
library(oligo)
library(sva)
library(pd.hugene.2.0.st)
library(matrixStats)
library(biomaRt)
library(WGCNA)
library(VennDiagram)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
library(ggplot2)

# Data ####
options(stringsAsFactors = FALSE)
Sys.setenv(R_THREADS = 1)
enableWGCNAThreads()
allowWGCNAThreads()

# Read in normalized counts #

exp_femdata <- read.csv("fcountsnu2.csv")
exp_maledata <- read.csv("mcountsnu2.csv")

getwd()

# Format normalized counts for WGCNA (input requires transposing) #

exp <- as.data.frame(t(exp_maledata[,-c(1)]))
names(exp) = exp_maledata$Gene_ID
rownames(exp) = names(exp_maledata)[-c(1)]
checkSets(exp)

exp2 <- list(maledata = list(data = as.data.frame(t(exp_maledata[,-c(1)]))))
names(exp2$maledata$data) = exp_maledata$Gene_ID
rownames(exp2$maledata$data) = names(exp_maledata)[-c(1)]
checkSets(exp2)

exp <- list(femdata = list(data = as.data.frame(t(exp_femdata[,-c(1)]))),
            maledata = list(data = as.data.frame(t(exp_maledata[,-c(1)]))))
names(exp$femdata$data) = exp_femdata$Gene_ID
rownames(exp$femdata$data) = names(exp_femdata)[-c(1)]
names(exp$maledata$data) = exp_maledata$Gene_ID
rownames(exp$maledata$data) = names(exp_maledata)[-c(1)]
checkSets(exp)

# Run WGCNA analysis #

mMods <- blockwiseModules(exp, checkMissingData = FALSE, maxBlockSize = 89, corType = "bicor",
                         maxPOutliers = 0.1, power = 5, networkType = "signed", 
                         checkPower = FALSE, TOMType = "signed", 
                         networkCalibration = "full quantile", saveTOMs = TRUE,
                         deepSplit = 4, mergeCutHeight = 0.1, verbose = 5)
table(mMods$colors) %>% sort(decreasing = TRUE)
module.dist <- as.data.frame(table(mMods$colors) %>% sort(decreasing = TRUE))
colnames(module.dist) <- c("Module", "Genes")
write.csv(module.dist,"module.distribution_phoapho.csv")

fMods <- blockwiseModules(exp$femdata$data, checkMissingData = FALSE, maxBlockSize = 50000, corType = "bicor",
                         maxPOutliers = 0.1, power = 6, networkType = "signed", 
                         checkPower = FALSE, minModuleSize = 100, TOMType = "signed", 
                         networkCalibration = "full quantile", saveTOMs = TRUE,
                         deepSplit = 4, mergeCutHeight = 0.1, verbose = 5)
table(fMods$colors) %>% sort(decreasing = TRUE)
female_module.dist <- as.data.frame(table(fMods$colors) %>% sort(decreasing = TRUE))
colnames(female_module.dist) <- c("Module", "Genes")

# Plot Merged Gene Dendrogram with Modules make sure to rename files #

pdf("males_dendogram.pdf", width = 10, height = 5)
sizeGrWindow(10, 5)
plotDendroAndColors(dendro = mMods$dendrograms[[1]], colors = mMods$colors, 
                    groupLabels = "Modules", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05, marAll = c(1, 5, 1, 0), main = "", cex.colorLabels = 1.3)
dev.off()

pdf("females_dendogram.pdf", width = 10, height = 5)
sizeGrWindow(10, 5)
plotDendroAndColors(dendro = fMods$dendrograms[[1]], colors = fMods$colors, 
                    groupLabels = "Modules", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05, marAll = c(1, 5, 1, 0), main = "", cex.colorLabels = 1.3)
dev.off()

# Cluster Modules by MALE Eigengenes and Plot Dendrogram
METree <- (1 - bicor(mMods$MEs, maxPOutliers = 0.1)) %>% as.dist %>% 
  hclust(method = "average")
pdf("MALE Module Eigengene Dendrogram.pdf", height = 5, width = 10)
sizeGrWindow(height = 5, width = 10)
par(mar = c(0, 5, 1, 1))
plot(METree, main = "", xlab = "", sub = "", ylim = c(0, 1), cex = 0.6)
abline(h = 0.1, col = "red")
dev.off()
rm(METree)

METree <- (1 - bicor(fMods$MEs, maxPOutliers = 0.1)) %>% as.dist %>% 
  hclust(method = "average")
pdf("FEMALE Module Eigengene Dendrogramnu4.pdf", height = 5, width = 10)
sizeGrWindow(height = 5, width = 10)
par(mar = c(0, 5, 1, 1))
plot(METree, main = "", xlab = "", sub = "", ylim = c(0, 1), cex = 0.6)
abline(h = 0.1, col = "red")
dev.off()
rm(METree)

# Compare Eigengene Networks Between FEMALE and MALE
MEs <- orderMEs(Mods$MEs)
pdf(file = "Male comparison WGCNA Eigengene Networksnu3.pdf", width = 8, height = 7)
sizeGrWindow(width = 8, height = 7)
par(cex = 0.8)
plotEigengeneNetworks(consensusMEs, setLabels = c("Male Cortex"), 
                      plotDendrograms = FALSE, marHeatmap = c(3, 3, 2, 1), zlimPreservation = c(0.5, 1), 
                      xLabelsAngle = 90)
dev.off()



# Calculate Module Membership ####
#moduleMembership <- mapply(bicorAndPvalue, exp, MEs, 
#                               MoreArgs = list(alternative = "two.sided", use = "pairwise.complete.obs", 
#                                               maxPOutliers = 0.1))
moduleMembership <- bicorAndPvalue(exp, MEs, use = "pairwise.complete.obs", alternative = "two.sided")


MM_male <- as.data.frame(moduleMembership$bicor)
colnames(MM_male) <- gsub(pattern = "ME", replacement = "", x = colnames(MM_male), fixed = TRUE)
MM_male$Probe <- rownames(MM_male)
MM_male$Module <- Mods$colors
write.table(MM_male, "Consensus Modules MALES Probe Module Membershipnu4.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# Get Module Hub Probes and Genes ####
hubProbes_male <- sapply(colnames(MM_male)[!colnames(MM_male) %in% c("Probe", "Module")], function(x){
  temp <- MM_male[MM_male$Module == x,]
  temp$Probe[temp[, x] == max(temp[, x])] %>% as.character %>% unique %>% sort})

ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

hubGenes_male <- lapply(hubProbes_male, function(x){
  getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", values = x, mart = ensembl, 
        verbose = TRUE) %>% unlist %>% as.character %>% unique %>% sort %>% paste(collapse = ", ")}) %>% unlist

# Combine Covariates ####
cov_male <- read.csv("MaleTraits.csv")
cov_male <- cov_male[,c("Mice", "Timepoint", "Light")]
colnames(cov_male)[colnames(cov_male) == "Mice"] <- "sampleID"

rownames(cov_male) <- cov_male$sampleID
cov_male <- cov_male[,c("Timepoint", "Light")]
cov_male <- as.matrix(cov_male)
table(rownames(exp$maledata$data) == rownames(cov_male)) # All TRUE
pheno <- list(male = list(data = cov_male))

# Get Meta-Analysis Correlations ####
moduleMembership <- read.delim("Consensus Modules MALES Probe Module Membershipnu4.txt", sep = "\t",
                               header = TRUE, stringsAsFactors = FALSE)
Mods2 <- moduleMembership$Module
MEs_male <- moduleEigengenes(t(exp_maledata[,-c(1)]), colors = Mods2)$eigengenes
rownames(MEs_male) <- rownames(t(exp_maledata[,-c(1)]))
consensusMEs <- list(male = list(data = MEs_male))
consensusMEs <- orderMEs(consensusMEs)

MEMAs <- list()
for (t in 1:20){
  MEMAs[[t]] = standardScreeningNumericTrait(MEs_male, subset(pheno$male$data, colIndex = t), 
                                             corFnc = bicor,
                                             corOptions = list(maxPOutliers = 0.1, use = "pairwise.complete.obs"), 
                                             qValues = TRUE,)
}

zscores <- sapply(MEMAs, function(x) x[["Z.RootDoFWeights"]])
rownames(zscores) <- colnames(consensusMEs$female$data)
colnames(zscores) <- colnames(pheno$female$data)
pvalues <- sapply(MEMAs, function(x) x[["p.RootDoFWeights"]])
dimnames(pvalues) <- dimnames(zscores)
qvalues <- sapply(MEMAs, function(x) x[["q.RootDoFWeights"]])
dimnames(qvalues) <- dimnames(zscores)

# Plot Correlations ####
# Plot All Correlations for Meta Analysis (Z-scores)
star <- apply(qvalues, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Consensus Modules Meta Covariate Correlation Plotnu3.pdf", width = 11, height = 15)
sizeGrWindow(width = 11, height = 15)
par(mar = c(9, 8, 1, 2))
labeledHeatmap(Matrix = zscores, xLabels = colnames(zscores), yLabels = rownames(zscores), 
               ySymbols = gsub("ME", "", rownames(zscores)), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = FALSE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-5, 5), main = "", cex.lab.y = 1)
dev.off()

# Plot Significant Correlations for Meta Analysis (Z-scores)
sigRows <- rowSums2(qvalues < 0.05) >= 1
sigCols <- colSums2(qvalues < 0.05) >= 1
zscores_sub <- zscores[sigRows, sigCols]
qvalues_sub <- qvalues[sigRows, sigCols]
colnames(zscores_sub) <- c("Timepoint", "Light")
star <- apply(qvalues_sub, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Consensus Modules Meta Covariate Correlation Plot Significant Onlynu3.pdf", width = 8, height = 11)
sizeGrWindow(width = 8, height = 11)
par(mar = c(6, 8, 1, 1))
labeledHeatmap(Matrix = zscores_sub, xLabels = colnames(zscores_sub), yLabels = rownames(zscores_sub), 
               ySymbols = gsub("ME", "", rownames(zscores_sub)), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = FALSE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-5, 5), main = "", cex.lab.y = 1)
dev.off()

# Plot All Correlations for FEMALE (Z-scores)
zscores <- sapply(MEMAs, function(x) x[["Z.female"]])
rownames(zscores) <- colnames(consensusMEs$female$data)
colnames(zscores) <- colnames(pheno$female$data)
pvalues <- sapply(MEMAs, function(x) x[[which(names(x) %in% c("pvalueStudent.female", "pvalueStudent.1.vs.2.female"))]])
dimnames(pvalues) <- dimnames(zscores)
qvalues <- sapply(MEMAs, function(x) x[[which(names(x) %in% c("qvalueStudent.female", "q.Student.female"))]])
dimnames(qvalues) <- dimnames(zscores)
star <- apply(qvalues, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Consensus Modules FEMALES Covariate Correlation Plot zscoresnu3.pdf", width = 11, height = 15)
sizeGrWindow(width = 11, height = 15)
par(mar = c(9, 8, 1, 2))
labeledHeatmap(Matrix = zscores, xLabels = colnames(zscores), yLabels = rownames(zscores), 
               ySymbols = gsub("ME", "", rownames(zscores)), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = FALSE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-5, 5), main = "", cex.lab.y = 1)
dev.off()

# Plot All Correlations for MALE (Z-scores)
zscores <- sapply(MEMAs, function(x) x[["Z.male"]])
rownames(zscores) <- colnames(consensusMEs$male$data)
colnames(zscores) <- colnames(pheno$male$data)
pvalues <- sapply(MEMAs, function(x) x[[which(names(x) %in% c("pvalueStudent.male", "pvalueStudent.1.vs.2.male"))]])
dimnames(pvalues) <- dimnames(zscores)
qvalues <- sapply(MEMAs, function(x) x[[which(names(x) %in% c("qvalueStudent.male", "q.Student.male"))]])
dimnames(qvalues) <- dimnames(zscores)
star <- apply(qvalues, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Consensus Modules MALE Covariate Correlation Plot zscoresnu3.pdf", width = 11, height = 15)
sizeGrWindow(width = 11, height = 15)
par(mar = c(9, 8, 1, 2))
labeledHeatmap(Matrix = zscores, xLabels = colnames(zscores), yLabels = rownames(zscores), 
               ySymbols = gsub("ME", "", rownames(zscores)), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = FALSE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-5, 5), main = "", cex.lab.y = 1)
dev.off()

# Put Together Meta-Analysis Stats Table ####
n_Probes <- sapply(gsub(pattern = "ME", replacement = "", 
                        x = colnames(MM_female)[!colnames(MM_female) %in% c("Probe", "Module")], fixed = TRUE),
                   function(x) length(MM_female$Module[MM_female$Module == x]))
moduleStats_meta <- cbind(colnames(consensusMEs$female$data),
                          n_Probes,
                          sapply(hubProbes_female, paste, collapse = ", "),
                          sapply(hubProbes_male, paste, collapse = ", "),
                          hubGenes_female,
                          hubGenes_male,
                          sapply(MEMAs, function(x) x[, "Z.female"]),
                          sapply(MEMAs, function(x) x[, colnames(x)[colnames(x) %in% c("pvalueStudent.female", "pvalueStudent.1.vs.2.female")]]),
                          sapply(MEMAs, function(x) x[, colnames(x)[colnames(x) %in% c("qvalueStudent.female", "q.Student.female")]]),
                          sapply(MEMAs, function(x) x[, "Z.male"]),
                          sapply(MEMAs, function(x) x[, colnames(x)[colnames(x) %in% c("pvalueStudent.male", "pvalueStudent.1.vs.2.male")]]),
                          sapply(MEMAs, function(x) x[, colnames(x)[colnames(x) %in% c("qvalueStudent.male", "q.Student.male")]]),
                          sapply(MEMAs, function(x) x[, "Z.RootDoFWeights"]),
                          sapply(MEMAs, function(x) x[, "p.RootDoFWeights"]),
                          sapply(MEMAs, function(x) x[, "q.RootDoFWeights"]))
moduleStats_meta <- as.data.frame(moduleStats_meta)
colnames(moduleStats_meta) <- c("Module", "n_Probes", "hubProbes_female", "hubProbes_male", "hubGenes_female", "hubGenes_male",
                                paste(rep(c("Zscore_female", "pvalue_female", "qvalue_female", "Zscore_male", "pvalue_male",
                                            "qvalue_male", "Zscore_Meta", "pvalue_Meta", "qvalue_Meta"), 
                                          each = length(colnames(pheno$female$data))), 
                                      rep(colnames(pheno$female$data), 5), sep = "_"))
write.table(moduleStats_meta, "Consensus Modules Meta Covariate Correlation Stats with Hub Genesnu3.txt", sep = "\t", quote = FALSE)



