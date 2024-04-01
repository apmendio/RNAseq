# Packages ####
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/WGCNA_Final")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/individualWGCNA/optimization_folder/test/comparative")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/individualWGCNA/optimization_folder/test/comparative/test-mb")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/clamsrw/rnaseq")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/clamsrw/rnaseq/7")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/clamsrw/rnaseq/7/combined")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/clamsrw/rnaseq/combined_analysis")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/clamsrw/rnaseq/combined_analysis/separated")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/clamsrw/rnaseq/combined_analysis/separated/test")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/clamsrw/rnaseq/separate_normalization")

sapply(c("tidyverse", "oligo", "sva","pd.hugene.2.0.st", "matrixStats", "biomaRt", "WGCNA", "VennDiagram"), require, character.only = TRUE)
getwd()
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

exp_femdata <- read.csv("baboon_allgenes_counts.csv")
exp_femdata2 <- read.csv("baboon_allgenes_counts2.csv")
exp_maledata <- read.csv("mcountsnu2.csv")
musdata <- exp_maledata
papdata <- exp_femdata
rownames(musdata) <- exp_maledata$Gene_ID
musdata$Gene_ID <- musgenes

#clams-rw

exp_femdata <- read.csv("clamsrw.csv")
exp_maledata <- read.csv("mcountsnu2.csv")

#male<-as.data.frame(exp_maledata)
#rownames(male$Gene_ID) = names(male)[-c(1)]
#names(male$Gene_ID) = male$Gene_ID
#female<-data.matrix(exp_femdata)
# Check directories #
musens <- exp_maledata$Gene_ID
musgenes2 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "refseq_mrna", "refseq_mrna_predicted"), filters = "ensembl_gene_id", values = exp_maledata$Gene_ID, mart = ensembl)
musgenes2 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters = "ensembl_gene_id", values = exp_maledata$Gene_ID, mart = ensembl)
colnames(musgenes2)[1] <- "Gene_ID"
test <- merge(musgenes2, musdata, by = "ensembl_gene_id")
colnames(test)[2]
colnames(test)[2] <- "Gene_ID"
colnames(musgenes2)[2] <- "Gene_ID"
colnames(musgenes2)[1] <- "ensemble_gene_id"
musdata <- test
colnames(musdata)[1] <- "ensembl_gene_id"
colnames(musdata)[2] <- "Gene_ID"
colnames(papdata)[1]
rownames(musdata) <- musdata$Gene_ID
mousedata = musdata
mousedata <- musdata[-c(1)]
papdata <- data.frame(papdata)
rn <- rownames(papdata) 
papdata$Gene_ID <- rn

library(dplyr)
test4 <- data.frame(tolower(test$Gene_ID))
test$test <- test4$tolower.test.Gene_ID.
test5 <- data.frame(tolower(papdata$Gene_ID))
papdata$test <- test5$tolower.papdata.Gene_ID.
test6 <- merge(test, papdata, by = "test")
test3 <- match(test$Gene_ID, papdata$Gene_ID)
test7<- data.frame(match(test$test, papdata$test))
test$Gene_ID

write.csv(test6, "combinedmb.csv")
getwd()

test_merge <- merge(exp_maledata, exp_femdata, by = "Gene_ID")
write.csv(test_merge, "combinedmb.csv")

exp_femdata2 <- read.csv("clams-rw_subset.csv")
exp_maledata2 <- read.csv("wtmalemouse.csv")
exp_femdata2 <- read.csv("fcountsnu2.csv")

cf <- read.csv("clams-cf_normalized.csv")
cm <- read.csv("clams-cm_normalized.csv")
rf <- read.csv("rw-rf_normalized.csv")
cfcm <- merge(cf, cm, by = "EnsemblID")
rm <- read.csv("rw-rm_normalized.csv")
wtfm <- read.csv("wtfm.csv")
merge <- merge(rm, exp_femdata2, by = "Gene_ID")
merge <- merge(cfcm, exp_maledata2, by = "EnsemblID")
cfcmrm <- read.csv("cfcmrm.csv")
write.csv(merge, "rf_traits.csv")

##2.20.24

cf <- read.csv("clams-cf_normalized.csv")
cm <- read.csv("clams-cm_normalized.csv")
rf <- read.csv("rw-rf_normalized.csv")
rm <- read.csv("rw-rm_normalized.csv")
wtfm <- read.csv("wtfm.csv")
merge <- merge(rm, exp_femdata2, by = "Gene_ID")
merge <- merge(cfcm, exp_maledata2, by = "EnsemblID")
cfcmrm <- read.csv("cfcmrm.csv")
write.csv(merge, "rf_traits.csv")

###
merge <- merge(exp_maledata2, exp_femdata2, by = "Gene_ID")
merge <- merge(wtfm, exp_femdata2, by = "Gene_ID")
write.csv(merge, "wtfm.csv")

exp_femdata2 <- read.csv("clams-rw_subset.csv")
exp_maledata2 <- read.csv("wtfm.csv")

exp_femdata2 <- read.csv("clams-rw_subset.csv")
exp_maledata2 <- read.csv("wtmalemouse.csv")
exp_femdata2 <- read.csv("fcountsnu2.csv")

merge <- merge(exp_maledata2, exp_femdata2, by = "Gene_ID")
write.csv(merge, "wtfm.csv")

exp_femdata2 <- read.csv("clams-rw_subset.csv")
exp_maledata2 <- read.csv("wtfm.csv")

exp_femdata2 <- read.csv("cfcmrm.csv")
exp_maledata2 <- read.csv("wtfm.csv")

##Separated analysis

exp_femdata2 <- read.csv("clams_subset.csv")
exp_femdata2 <- read.csv("rw_subset.csv")
exp_maledata2 <- read.csv("wtfm.csv")

#Separate Normalization test
exp_femdata2 <- read.csv("cfcm.csv")
exp_maledata2 <- read.csv("wtfm.csv")

exp_femdata2 <- read.csv("cfcmrm.csv")
exp_maledata2 <- read.csv("wtfm.csv")
# Format normalized counts for WGCNA (input requires transposing) #

#head(exp_femdata[,-c(1)])
exp <- list(femdata = list(data = as.data.frame(t(exp_femdata[,-c(1)]))),
            maledata = list(data = as.data.frame(t(exp_maledata[,-c(1)]))))
names(exp$femdata$data) = exp_femdata$Gene_ID
rownames(exp$femdata$data) = names(exp_femdata)[-c(1)]
names(exp$maledata$data) = exp_maledata$Gene_ID
rownames(exp$maledata$data) = names(exp_maledata)[-c(1)]
checkSets(exp)

#clams-rw
#head(exp_femdata[,-c(1)])
exp <- list(femdata = list(data = as.data.frame(t(exp_femdata2[,-c(1)]))),
            maledata = list(data = as.data.frame(t(exp_maledata2[,-c(1)]))))
names(exp$femdata$data) = exp_femdata2$Gene_ID
rownames(exp$femdata$data) = names(exp_femdata2)[-c(1)]
names(exp$maledata$data) = exp_maledata2$Gene_ID
rownames(exp$maledata$data) = names(exp_maledata2)[-c(1)]
checkSets(exp)
exp

exp <- list(femdata = list(data = as.data.frame(t(exp_femdata2[,-c(1)]))),
            maledata = list(data = as.data.frame(t(exp_maledata2[,-c(1)]))))
names(exp$femdata$data) = exp_femdata2$EnsemblID
rownames(exp$femdata$data) = names(exp_femdata2)[-c(1)]
names(exp$maledata$data) = exp_maledata2$EnsemblID
rownames(exp$maledata$data) = names(exp_maledata2)[-c(1)]
checkSets(exp)
exp

###2.20.24

wtfm_norm <- read.csv("wtfm_normalized.csv")
cf_norm <- read.csv("cf_normalized.csv")
cm_norm <- read.csv("cm_normalized.csv")
rf_norm <- read.csv("rf_normalized.csv")
rm_norm <- read.csv("rm_normalized.csv")

row.names(wtfm_norm) <- wtfm_norm$Gene_ID
row.names(cf_norm) <- cf_norm$Gene_ID
row.names(cm_norm) <- cm_norm$Gene_ID
row.names(rf_norm) <- rf_norm$Gene_ID
row.names(rm_norm) <- rm_norm$Gene_ID

multiExpr <- list(wtfdata = list(data = as.data.frame(t(wtfm_norm[,-c(1)]))),
                  wtmdata = list(data = as.data.frame(t(cf_norm[,-c(1)]))),
                  cfdata = list(data = as.data.frame(t(cf_norm[,-c(1)]))),
                  cmdata = list(data = as.data.frame(t(cm_norm[,-c(1)]))),
                  rfdata = list(data = as.data.frame(t(rf_norm[,-c(1)]))),
                  rmdata = list(data = as.data.frame(t(rm_norm[,-c(1)]))))

multiExpr <- list(wtfdata = list(data = as.data.frame(t(wtf_norm))),
                  wtmdata = list(data = as.data.frame(t(wtm_norm))),
                  cfdata = list(data = as.data.frame(t(cf_norm))),
                  cmdata = list(data = as.data.frame(t(cm_norm))),
                  rfdata = list(data = as.data.frame(t(rf_norm))),
                  rmdata = list(data = as.data.frame(t(rm_norm))))
checkSets(multiExpr)
exp <- list(cfdata = list(data = as.data.frame(t(cf_norm[,-c(1)]))),
                  cmdata = list(data = as.data.frame(t(cm_norm[,-c(1)]))))

exp2 <- list(rfdata = list(data = as.data.frame(t(rf_norm[,-c(1)]))),
                  rmdata = list(data = as.data.frame(t(rm_norm[,-c(1)]))))

# Run WGCNA analysis #

consensusMods <- blockwiseConsensusModules(exp, checkMissingData = FALSE, maxBlockSize = 14000, corType = "bicor",
                                           maxPOutliers = 0.1, power = 14, networkType = "signed", 
                                           checkPower = FALSE, minModuleSize = 30, TOMType = "signed", 
                                           networkCalibration = "full quantile", saveConsensusTOMs = TRUE,
                                           deepSplit = 4, mergeCutHeight = 0.1, verbose = 5)
table(consensusMods$colors) %>% sort(decreasing = TRUE)
module.dist <- as.data.frame(table(consensusMods$colors) %>% sort(decreasing = TRUE))
colnames(module.dist) <- c("Module", "Genes")
write.csv(module.dist,"module.distribution_consensusclams_signed.csv")

# Plot Merged Gene Dendrogram with Modules make sure to rename files #

pdf("consensus_wgcna_clams.pdf", width = 10, height = 5)
sizeGrWindow(10, 5)
plotDendroAndColors(dendro = consensusMods$dendrograms[[1]], colors = consensusMods$colors, 
                    groupLabels = "Modules", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05, marAll = c(1, 5, 1, 0), main = "", cex.colorLabels = 1.3)
dev.off()

# Cluster Modules by FEMALE Eigengenes and Plot Dendrogram
METree <- (1 - bicor(consensusMods$multiMEs$cfdata$data, maxPOutliers = 0.1)) %>% as.dist %>% 
  hclust(method = "average")
pdf("CLAMS FEMALE Module Eigengene Dendrogramnu.pdf", height = 5, width = 10)
sizeGrWindow(height = 5, width = 10)
par(mar = c(0, 5, 1, 1))
plot(METree, main = "", xlab = "", sub = "", ylim = c(0, 1), cex = 0.6)
abline(h = 0.1, col = "red")
dev.off()

# Cluster Modules by MALE Eigengenes and Plot Dendrogram
METree <- (1 - bicor(consensusMods$multiMEs$cmdata$data, maxPOutliers = 0.1)) %>% as.dist %>% 
  hclust(method = "average")
pdf("CLAMS MALE Module Eigengene Dendrogramnu.pdf", height = 5, width = 10)
sizeGrWindow(height = 5, width = 10)
par(mar = c(0, 5, 1, 1))
plot(METree, main = "", xlab = "", sub = "", ylim = c(0, 1), cex = 0.6)
abline(h = 0.1, col = "red")
dev.off()
rm(METree)

# Compare Eigengene Networks Between FEMALE and MALE
consensusMEs <- consensusOrderMEs(consensusMods$multiMEs)
pdf(file = "CLAMS Consensus WGCNA Eigengene Networks.pdf", width = 8, height = 7)
sizeGrWindow(width = 8, height = 7)
par(cex = 0.8)
plotEigengeneNetworks(consensusMEs, setLabels = c("CLAMS - Female Cortex", "CLAMS - Male Cortex"), 
                      plotDendrograms = FALSE, marHeatmap = c(3, 3, 2, 1), zlimPreservation = c(0.5, 1), 
                      xLabelsAngle = 90)
dev.off()

# Calculate Module Membership ####
moduleMembership2 <- mtd.mapply(corAndPvalue, exp, consensusMEs, 
                               MoreArgs = list(alternative = "two.sided", use = "pairwise.complete.obs")) 
                                               maxPOutliers = 0.1))
MM_female <- as.data.frame(moduleMembership$cfdata$data$bicor)
colnames(MM_female) <- gsub(pattern = "ME", replacement = "", x = colnames(MM_female), fixed = TRUE)
MM_female$Probe <- rownames(MM_female)
MM_female$Module <- consensusMods$colors
write.table(MM_female, "Consensus Modules FEMALES CLAMS Module Membership.txt", sep = "\t", quote = FALSE, row.names = FALSE)

MM_male <- as.data.frame(moduleMembership$cmdata$data$bicor)
colnames(MM_male) <- gsub(pattern = "ME", replacement = "", x = colnames(MM_male), fixed = TRUE)
MM_male$Probe <- rownames(MM_male)
MM_male$Module <- consensusMods$colors
write.table(MM_male, "Consensus Modules MALES CLAMS Module Membershipn.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Get Module Hub Probes and Genes ####
hubProbes_female <- sapply(colnames(MM_female)[!colnames(MM_female) %in% c("Probe", "Module")], function(x){
  temp <- MM_female[MM_female$Module == x,]
  temp$Probe[temp[, x] == max(temp[, x])] %>% as.character %>% unique %>% sort})
hubProbes_male <- sapply(colnames(MM_male)[!colnames(MM_male) %in% c("Probe", "Module")], function(x){
  temp <- MM_male[MM_male$Module == x,]
  temp$Probe[temp[, x] == max(temp[, x])] %>% as.character %>% unique %>% sort})
ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
hubGenes_female <- lapply(hubProbes_female, function(x){
  getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", values = x, mart = ensembl, 
        verbose = TRUE) %>% unlist %>% as.character %>% unique %>% sort %>% paste(collapse = ", ")}) %>% unlist
hubGenes_female
write.csv(hubGenes_female, "hubGenes_CLAMS-RW9.csv")
hubGenes_male <- lapply(hubProbes_male, function(x){
  getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", values = x, mart = ensembl, 
        verbose = TRUE) %>% unlist %>% as.character %>% unique %>% sort %>% paste(collapse = ", ")}) %>% unlist
hubGenes_male
write.csv(hubGenes_male, "hubGenes_WT9.csv")

### CLAMS-RW ###
# Combine Covariates ####
cov_female <- read.csv("cf_traits.csv")
cov_male <- read.csv("cm_traits.csv")
cov_female <- cov_female[,c("SampleID", "Timepoint", "Light")]
colnames(cov_female)[colnames(cov_female) == "SampleID"] <- "SampleID"
cov_male <- cov_male[,c("Mice", "Timepoint", "Light")]
colnames(cov_male)[colnames(cov_male) == "Mice"] <- "SampleID"
table(colnames(cov_female) == colnames(cov_male)) # All TRUE

rownames(cov_female) <- cov_female$SampleID
#cov_female <- cov_female[,c("Timepoint", "Light")]
cov_female <- as.matrix(cov_female)
cov_female <- cov_female[,-c(1,4,5,6,7,8)]
table(rownames(exp$cf$data) == rownames(cov_female)) # All TRUE

rownames(cov_male) <- cov_male$SampleID
#cov_male <- cov_male[,c("Timepoint", "Light")]
#cov_male <- as.matrix(cov_male)
cov_male <- cov_male[,-c(1,4,5,6,7,8)]
table(rownames(exp$cmdata$data) == rownames(cov_male)) # All TRUE
pheno <- list(female = list(data = cov_female),
              male = list(data = cov_male))

### Module Trait Correlation

moduleTraitCorP = cor(consensusMEs$cfdata$data, pheno$female$data, method = "pearson", use = "p")
moduleTraitPvalueP = corPvalueFisher(moduleTraitCorP, 32, twoSided = TRUE)

moduleTraitCorS = cor(consensusMEs$cfdata$data, pheno$female$data, method = "spearman", use = "p")
moduleTraitPvalueS = corPvalueFisher(moduleTraitCorS, 32, twoSided = TRUE)

# Get Meta-Analysis Correlations ####

MEMAs <- list()
for (t in 1:20){
  MEMAs[[t]] = metaAnalysis(consensusMEs, mtd.subset(pheno, colIndex = t), useRankPvalue = FALSE, 
                            corFnc = bicor,
                            corOptions = list(maxPOutliers = 0.1, use = "pairwise.complete.obs"), 
                            getQvalues = TRUE)
}

zscores <- sapply(MEMAs, function(x) x[["Z.RootDoFWeights"]])
rownames(zscores) <- colnames(consensusMEs$cfdata$data)
colnames(zscores) <- colnames(pheno$female$data)
pvalues <- sapply(MEMAs, function(x) x[["p.RootDoFWeights"]])
dimnames(pvalues) <- dimnames(zscores)
qvalues <- sapply(MEMAs, function(x) x[["q.RootDoFWeights"]])
dimnames(qvalues) <- dimnames(zscores)

# Plot Correlations ####
# Plot All Correlations for Meta Analysis (Z-scores)
star <- apply(qvalues, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Consensus Modules Correlation Plot CLAMS Female Pearson.pdf", width = 11, height = 15)
sizeGrWindow(width = 11, height = 15)
par(mar = c(9, 8, 1, 2))
labeledHeatmap(Matrix = zscores, xLabels = colnames(zscores), yLabels = rownames(zscores), 
               ySymbols = gsub("ME", "", rownames(zscores)), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = FALSE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-5, 5), main = "", cex.lab.y = 1, plotLegend = TRUE, legendLabel = "Z-Score")
dev.off()

# Plot Correlations for individual datasets

star <- apply(moduleTraitPvalueP, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Consensus Modules CLAMS Female Pearson Correlation Plot.pdf", width = 11, height = 15)
sizeGrWindow(width = 11, height = 15)
par(mar = c(9, 8, 1, 2))
labeledHeatmap(Matrix = moduleTraitCorP, xLabels = colnames(zscores), yLabels = rownames(zscores), 
               ySymbols = gsub("ME", "", rownames(zscores)), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = TRUE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-1, 1), main = "Module Trait Correlations Pearson Correlation", cex.lab.y = 1, plotLegend = TRUE, legendLabel = "Pearson Correlation Coefficient")
dev.off()

star <- apply(moduleTraitPvalueS, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Consensus Modules CLAMS Female Spearman Correlation Plot.pdf", width = 11, height = 15)
sizeGrWindow(width = 11, height = 15)
par(mar = c(9, 8, 1, 2))
labeledHeatmap(Matrix = moduleTraitCorS, xLabels = colnames(zscores), yLabels = rownames(zscores), 
               ySymbols = gsub("ME", "", rownames(zscores)), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = TRUE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-1, 1), main = "Module Trait Correlations Spearman Correlation", cex.lab.y = 1, plotLegend = TRUE, legendLabel = "Spearman Correlation Coefficient")
dev.off()

library(WGCNA)
star <- apply(moduleTraitPvaluecf, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
sizeGrWindow(width = 11, height = 15)
par(mar = c(9, 8, 1, 2))
labeledHeatmap(Matrix = moduleTraitCorcf, xLabels = colnames(moduleTraitCorcf), yLabels = rownames(moduleTraitCorcf), 
               ySymbols = gsub("ME", "", rownames(moduleTraitCorcf)), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = TRUE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-1, 1), main = "Module Trait Correlations Pearson Correlation", cex.lab.y = 1, plotLegend = TRUE, legendLabel = "Pearson Correlation Coefficient")
dev.off()

star <- apply(moduleTraitPvalueS2, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Consensus Modules Meta Covariate Correlation Plot RM-S.pdf", width = 11, height = 15)
sizeGrWindow(width = 11, height = 15)
par(mar = c(9, 8, 1, 2))
labeledHeatmap(Matrix = moduleTraitCorS, xLabels = colnames(zscores), yLabels = rownames(zscores), 
               ySymbols = gsub("ME", "", rownames(zscores)), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = TRUE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-1, 1), main = "Module Trait Correlations Spearman Correlation", cex.lab.y = 1, plotLegend = TRUE, legendLabel = "Spearman Correlation Coefficient")
dev.off()

star <- apply(moduleTraitPvalueK2, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Consensus Modules Meta Covariate Correlation Plot RW-K2-2s3.pdf", width = 11, height = 15)
sizeGrWindow(width = 11, height = 15)
par(mar = c(9, 8, 1, 2))
labeledHeatmap(Matrix = moduleTraitCorK, xLabels = colnames(zscores), yLabels = rownames(zscores), 
               ySymbols = gsub("ME", "", rownames(zscores)), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = TRUE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-1, 1), main = "Module Trait Correlations Kendall Correlation", cex.lab.y = 1, plotLegend = TRUE, legendLabel = "Kendall Correlation Coefficient")
dev.off()
# Plot Significant Correlations for Meta Analysis (Z-scores)
sigRows <- rowSums2(qvalues < 0.05) >= 1
sigCols <- colSums2(qvalues < 0.05) >= 1
zscores_sub <- zscores[sigRows, sigCols]
qvalues_sub <- qvalues[sigRows, sigCols]
colnames(zscores_sub) <- colnames(pheno$female$data)
star <- apply(qvalues_sub, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Consensus Modules Meta Covariate Correlation Plot Significant Only.pdf", width = 8, height = 11)
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
dimnames(moduleTraitPvalue) <- dimnames(zscores)
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

# Always change the color of module for Module Expression ####
# change the dataset between male and female to acquire the correct probes
honeydew_probes <- colnames(exp$femaledata$data)[consensusMods$colors == "honeydew"]
ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
honeydew_genes <- getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", 
                        values = honeydew_probes, mart = ensembl) %>% unlist %>% as.character %>% unique %>% sort
table(honeydew_genes)
write.table(honeydew_genes, "malehoneydewModule_geneslist.txt", sep = "\t", quote = FALSE)

##################################################Module gene list ####################################################
honeydew_probes_male <- colnames(exp$maledata$data)[consensusMods$colors == "orangered4"]
ensembl_male <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
honeydew_genes_male <- getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", 
                             values = honeydew_probes_male, mart = ensembl_male) %>% unlist %>% as.character %>% unique %>% sort
table(honeydew_genes_male)

honeydew_probes_female <- colnames(exp$femdata$data)[consensusMods$colors == "honeydew"]
ensembl_female <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
honeydew_genes_female <- getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", 
                               values = honeydew_probes_female, mart = ensembl) %>% unlist %>% as.character %>% unique %>% sort
table(honeydew_genes_female, honeydew_genes_male)
#######################################################################################################################
table(honeydew_genes)
# skymaroon1 module eigengene plots
maroon_ME_female <- consensusMEs$female$data$MEmaroon
maroon_ME_male <- consensusMEs$male$data$MEmaroon
maroon_ME <- as.data.frame(cbind(c(rownames(pheno$female$data), rownames(pheno$male$data)),
                                 c(pheno$female$data[, "Timepoint"], pheno$male$data[, "Timepoint"]),
                                 c(pheno$female$data[, "Light"], pheno$male$data[, "Light"]),
                                 c(maroon_ME_female, maroon_ME_male)))
colnames(maroon_ME) <- c("sampleID", "Timepoint", "Light", "ME")
maroon_ME$Timepoint[maroon_ME$Timepoint == 0] <- "ZT0"
maroon_ME$Timepoint[maroon_ME$Timepoint == 3] <- "ZT3"
maroon_ME$Timepoint[maroon_ME$Timepoint == 6] <- "ZT6"
maroon_ME$Timepoint[maroon_ME$Timepoint == 9] <- "ZT9"
maroon_ME$Timepoint[maroon_ME$Timepoint == 12] <- "ZT12"
maroon_ME$Timepoint[maroon_ME$Timepoint == 15] <- "ZT15"
maroon_ME$Timepoint[maroon_ME$Timepoint == 18] <- "ZT18"
maroon_ME$Timepoint[maroon_ME$Timepoint == 21] <- "ZT21"
maroon_ME$Light[maroon_ME$Light == 1] <- "ON"
maroon_ME$Light[maroon_ME$Light == 2] <- "OFF"
maroon_ME$Study <- c(rep("FEMALE", length(maroon_ME_female)), rep("MALE", length(maroon_ME_male)))
maroon_ME$Timepoint <- factor(maroon_ME$Timepoint, levels = c("ZT0", "ZT3", "ZT6", "ZT9", "ZT12", "ZT15", "ZT18", "ZT21"))
maroon_ME$Light <- factor(maroon_ME$Light, levels = c("ON", "OFF"))
maroon_ME$Study <- factor(maroon_ME$Study, levels = c("FEMALE", "MALE"))
maroon_ME$ME <- as.numeric(maroon_ME$ME)
maroon_ME <- maroon_ME[, c("sampleID", "Study", "Timepoint", "Light", "ME")]
write.csv(maroon_ME, "maroon_ME_traitsnu.csv")
maroon_ME <- read.csv("maroon_ME_traitsnu.csv")
rownames(maroon_ME) <- maroon_ME[,c(1)]
maroon_ME <- maroon_ME[,-c(1)]

gg <- ggplot(data = maroon_ME)
gg <- gg +
  geom_boxplot(aes(x = Light, y = ME, fill = Timepoint), size = 0.8, outlier.size = 0.8) +
  theme_bw(base_size = 24) +
  theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
        legend.key = element_blank(), panel.grid.minor = element_blank(),
        legend.position = c(0.87, 0.9), legend.background = element_blank(), 
        legend.key.size = unit(0.8, "cm"), strip.text.x = element_text(size = 22), 
        axis.ticks = element_line(size = 1.25), legend.title = element_text(size = 22),
        strip.background = element_blank(), legend.direction = "vertical", panel.spacing.y = unit(0, "lines"), 
        plot.margin = unit(c(0, 1, 1, 0.5), "lines"), axis.title.x = element_blank(), 
        axis.text = element_text(size = 22, color = "black")) +
  ylab("maroon Module Eigengene") +
  scale_fill_manual(name = "Timepoint", values = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  scale_color_manual(name = "Timepoint", values = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  coord_cartesian(ylim = c(-1, 1)) +
  facet_wrap(facets = vars(Study), nrow = 2)
ggsave(filename = "maroon Module Eigengene Boxplotnu.png", plot = gg, dpi = 600, units = "in")



