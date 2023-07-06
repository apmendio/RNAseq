setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/clamsrw/rnaseq/separate_normalization")
# module preservation
#mods = list();
# Sof thresholding powers for network definition
power = c(3, 2);
#collectGarbage();
#labels = list();
#nn = if (interactive()) nSets else 1;

exp_wtdata <- read.csv("wtfm.csv")
exp_cfdata <- read.csv("cf.csv")
exp_cmdata <- read.csv("cm.csv")
exp_rmdata <- read.csv("rm.csv")

multiExpr <- list(wtdata = list(data = as.data.frame(t(exp_wtdata[,-c(1)]))),
            cfdata = list(data = as.data.frame(t(exp_cfdata[,-c(1)]))),
            cmdata = list(data = as.data.frame(t(exp_cmdata[,-c(1)]))),
            rmdata = list(data = as.data.frame(t(exp_rmdata[,-c(1)]))))
names(multiExpr$wtdata$data) = exp_wtdata$Gene_ID
rownames(multiExpr$wtdata$data) = names(exp_wtdata)[-c(1)]
names(multiExpr$cfdata$data) = exp_cfdata$Gene_ID
rownames(multiExpr$cfdata$data) = names(exp_cfdata)[-c(1)]
names(multiExpr$cmdata$data) = exp_cmdata$Gene_ID
rownames(multiExpr$cmdata$data) = names(exp_cmdata)[-c(1)]
names(multiExpr$rmdata$data) = exp_rmdata$Gene_ID
rownames(multiExpr$rmdata$data) = names(exp_rmdata)[-c(1)]
checkSets(multiExpr)

wtmods <- blockwiseModules(multiExpr$wtdata$data, checkMissingData = FALSE, maxBlockSize = 14000, corType = "bicor",
                               maxPOutliers = 0.1, power = 5, networkType="signed",
                               checkPower = FALSE, minModuleSize = 100, TOMType = "signed",
                               networkCalibration = "full quantile", saveConsensusTOMs = TRUE,
                               deepSplit=4, mergeCutHeight = 0.1, verbose=4); 
table(wtmods$colors) %>% sort(decreasing = TRUE)
module.dist <- as.data.frame(table(wtmods$colors) %>% sort(decreasing = TRUE))
colnames(module.dist) <- c("Module", "Genes")
write.csv(module.dist,"module.distribution_wt0.8_sft5.csv")

#clamsrwmods <- blockwiseModules(exp$femdata$data, checkMissingData = FALSE, maxBlockSize = 14000, corType = "bicor",
                           # maxPOutliers = 0.1, power = 9, networkType="signed",
                           # checkPower = FALSE, minModuleSize = 100, TOMType = "signed",
                           # networkCalibration = "full quantile", saveConsensusTOMs = TRUE,
                           # deepSplit=4, mergeCutHeight = 0.1, verbose=4); 
#table(clamsrwmods$colors) %>% sort(decreasing = TRUE)
#module.dist <- as.data.frame(table(clamsrwmods$colors) %>% sort(decreasing = TRUE))
#colnames(module.dist) <- c("Module", "Genes")
#write.csv(module.dist,"module.distribution_clams-rw0.8_sft9.csv")

#clamsmods <- blockwiseModules(exp2$femdata$data, checkMissingData = FALSE, maxBlockSize = 14000, corType = "bicor",
                              #  maxPOutliers = 0.1, power = 8, networkType="signed",
                              #  checkPower = FALSE, minModuleSize = 100, TOMType = "signed",
                              #  networkCalibration = "full quantile", saveConsensusTOMs = TRUE,
                              #  deepSplit=4, mergeCutHeight = 0.1, verbose=4); 
#table(clamsmods$colors) %>% sort(decreasing = TRUE)
#module.dist <- as.data.frame(table(clamsmods$colors) %>% sort(decreasing = TRUE))
#colnames(module.dist) <- c("Module", "Genes")
#write.csv(module.dist,"module.distribution_clams-cfcm0.8_sft8.csv")

rwmods <- blockwiseModules(multiExpr$rmdata$data, checkMissingData = FALSE, maxBlockSize = 14000, corType = "bicor",
                                maxPOutliers = 0.1, power = 5, networkType="signed",
                                checkPower = FALSE, minModuleSize = 100, TOMType = "signed",
                                networkCalibration = "full quantile", saveConsensusTOMs = TRUE,
                                deepSplit=4, mergeCutHeight = 0.1, verbose=4); 
table(rwmods$colors) %>% sort(decreasing = TRUE)
module.dist <- as.data.frame(table(rwmods$colors) %>% sort(decreasing = TRUE))
colnames(module.dist) <- c("Module", "Genes")
write.csv(module.dist,"module.distribution_rm0.8_sft5.csv")

cfmods <- blockwiseModules(multiExpr$cfdata$data, checkMissingData = FALSE, maxBlockSize = 14000, corType = "bicor",
                           maxPOutliers = 0.1, power = 3, networkType="signed",
                           checkPower = FALSE, minModuleSize = 100, TOMType = "signed",
                           networkCalibration = "full quantile", saveConsensusTOMs = TRUE,
                           deepSplit=4, mergeCutHeight = 0.1, verbose=4); 
table(cfmods$colors) %>% sort(decreasing = TRUE)
module.dist <- as.data.frame(table(cfmods$colors) %>% sort(decreasing = TRUE))
colnames(module.dist) <- c("Module", "Genes")
write.csv(module.dist,"module.distribution_cf0.8_sft3.csv")

cmmods <- blockwiseModules(multiExpr$cmdata$data, checkMissingData = FALSE, maxBlockSize = 14000, corType = "bicor",
                           maxPOutliers = 0.1, power = 7, networkType="signed",
                           checkPower = FALSE, minModuleSize = 100, TOMType = "signed",
                           networkCalibration = "full quantile", saveConsensusTOMs = TRUE,
                           deepSplit=4, mergeCutHeight = 0.1, verbose=4); 
table(cmmods$colors) %>% sort(decreasing = TRUE)
module.dist <- as.data.frame(table(cmmods$colors) %>% sort(decreasing = TRUE))
colnames(module.dist) <- c("Module", "Genes")
write.csv(module.dist,"module.distribution_cm0.8_sft7.csv")

# Generate MEs
wtMEs <- orderMEs(wtmods$MEs)
#clamsrwMEs <- orderMEs(clamsrwmods$MEs)
#clamsMEs <- orderMEs(clamsmods$MEs)
rwMEs <- orderMEs(rwmods$MEs)
cfMEs <- orderMEs(cfmods$MEs)
cmMEs <- orderMEs(cmmods$MEs)

# Calculate Module Membership ####

moduleMembershipwt <- bicorAndPvalue(multiExpr$wtdata$data, wtMEs, 
                              alternative = "two.sided", use = "pairwise.complete.obs", 
                                              maxPOutliers = 0.1)

#moduleMembershipcrw <- bicorAndPvalue(exp$femdata$data, clamsrwMEs, 
                                     #alternative = "two.sided", use = "pairwise.complete.obs", 
                                     #maxPOutliers = 0.1)

#moduleMembershipcfcm <- bicorAndPvalue(exp2$femdata$data, clamsMEs, 
                                      #alternative = "two.sided", use = "pairwise.complete.obs", 
                                      #maxPOutliers = 0.1)

moduleMembershiprw <- bicorAndPvalue(multiExpr$rmdata$data, rwMEs, 
                                      alternative = "two.sided", use = "pairwise.complete.obs", 
                                      maxPOutliers = 0.1)

moduleMembershipcf <- bicorAndPvalue(multiExpr$cfdata$data, cfMEs, 
                                      alternative = "two.sided", use = "pairwise.complete.obs", 
                                      maxPOutliers = 0.1)

moduleMembershipcm <- bicorAndPvalue(multiExpr$cmdata$data, cmMEs, 
                                      alternative = "two.sided", use = "pairwise.complete.obs", 
                                      maxPOutliers = 0.1)

ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host="https://useast.ensembl.org")

# Generate module membership files

MM_wt <- as.data.frame(moduleMembershipwt$bicor)
colnames(MM_wt) <- gsub(pattern = "ME", replacement = "", x = colnames(MM_wt), fixed = TRUE)
MM_wt$Probe <- rownames(MM_wt)
MM_wt$Module <- wtmods$colors

write.table(MM_wt, "WT Module Membership.txt", sep = "\t", quote = FALSE, row.names = FALSE)

MM_rw <- as.data.frame(moduleMembershiprw$bicor)
colnames(MM_rw) <- gsub(pattern = "ME", replacement = "", x = colnames(MM_rw), fixed = TRUE)
MM_rw$Probe <- rownames(MM_rw)
MM_rw$Module <- rwmods$colors

write.table(MM_rw, "RW Module Membership.txt", sep = "\t", quote = FALSE, row.names = FALSE)

MM_cf <- as.data.frame(moduleMembershipcf$bicor)
colnames(MM_cf) <- gsub(pattern = "ME", replacement = "", x = colnames(MM_cf), fixed = TRUE)
MM_cf$Probe <- rownames(MM_cf)
MM_cf$Module <- cfmods$colors

write.table(MM_cf, "CF Module Membership.txt", sep = "\t", quote = FALSE, row.names = FALSE)

MM_cm <- as.data.frame(moduleMembershipcm$bicor)
colnames(MM_cm) <- gsub(pattern = "ME", replacement = "", x = colnames(MM_cm), fixed = TRUE)
MM_cm$Probe <- rownames(MM_cm)
MM_cm$Module <- cmods$colors

write.table(MM_cm, "CM Module Membership.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Get Module Hub Probes and Genes ####
hubProbes_female <- sapply(colnames(MM_female)[!colnames(MM_female) %in% c("Probe", "Module")], function(x){
  temp <- MM_female[MM_female$Module == x,]
  temp$Probe[temp[, x] == max(temp[, x])] %>% as.character %>% unique %>% sort})
hubProbes_male <- sapply(colnames(MM_male)[!colnames(MM_male) %in% c("Probe", "Module")], function(x){
  temp <- MM_male[MM_male$Module == x,]
  temp$Probe[temp[, x] == max(temp[, x])] %>% as.character %>% unique %>% sort})
ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host="https://useast.ensembl.org")
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
