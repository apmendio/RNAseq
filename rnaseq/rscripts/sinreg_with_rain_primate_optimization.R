#### Using sine regression to assess associations between WGCNA gene modules and time ####

library(WGCNA)
BiocManager::install("rain")
library(rain)
library(ShellChron)
library(ggplot2)
library(openxlsx)
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Mm.eg.db")
library(AnnotationDbi)
library(org.Mm.eg.db)
library(RColorBrewer)
library(data.table)
getwd()
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/individualWGCNA")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/individualWGCNA/optimization_folder/test")

MEnames = colnames(MEs_male)
no.MEs = length(MEnames)

traits = pheno$male$data

alldata = merge(traits, MEs_male, by="row.names")

alldata = alldata[
  with(alldata, order(alldata$Timepoint)),
]

head(alldata)
write.csv(alldata, "wt-mf_clams.ev.csv")

fitlist = as.list(1:no.MEs)
names(fitlist) <- MEnames

p.values = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(p.values) = MEnames
colnames(p.values) = "p.value"

R2 = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(R2) = MEnames
colnames(R2) = "R2"

period = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(period) = MEnames
colnames(period) = "period"

peak = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(peak) = MEnames
colnames(peak) = "peak"

amplitude = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(amplitude) = MEnames
colnames(amplitude) = "amplitude"

plotpoints = data.frame(matrix(ncol=no.MEs, nrow=nrow(alldata)))
rownames(plotpoints) = alldata$Row.names
colnames(plotpoints) = MEnames

for(i in MEnames){
  
  # print status
  print(paste("Running entity:", i, "which is", which(MEnames==i), "out of", no.MEs))
  
  #create temporary data matrix and model formula
  
  x = alldata$Timepoint
  y = alldata[,i]
  fitlist[[i]] <- sinreg(x, y, plot=FALSE)
  p.values[i,1] <- fitlist[[i]][[1]][6]
  R2[i,1] <- fitlist[[i]][[1]][5]
  period[i,1] <- fitlist[[i]][[1]][3]
  peak[i,1] <- fitlist[[i]][[1]][4]
  amplitude[i,1] <- fitlist[[i]][[1]][2]
  plotpoints[,i] <- fitlist[[i]][[2]]
}

results = cbind(p.values, R2, period, peak, amplitude)

results$FDR = p.adjust(results$p.value, method="fdr")
plotpoints$Timepoint = alldata$Timepoint

for(i in MEnames) {
  y = alldata[,i]
  y1 = plotpoints[,i]
  ggplot(data=alldata, aes(x=Timepoint, y=y)) +
    geom_point() +
    geom_line(data=plotpoints, aes(x=Timepoint, y=y1, color="red")) +
    scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21)) +
    ggtitle(paste(i)) +
    ylab("Module EigenValue") +
    xlab("Zeitgeber Time (ZT)") +
    theme_classic() +
    theme(legend.position="none")
  ggsave(paste("Sine_regression_plot_", i, "_WTfm_clams.pdf"))
}

# order results by best fit (highest R2 value) and save results #
results.ordered = results[
  with(results, order(results$R2, decreasing=TRUE)),
]

results.ordered$Module = rownames(results.ordered)

openxlsx::write.xlsx(results.ordered, file="Sine_regresion_results_WTfm_clams.xlsx")

# Checking if hub genes within modules cycle in a similar manner to EigenValues #

exp_maledata.1 = as.data.frame(t(exp_maledata2[,-1]))
colnames(exp_maledata.1) = exp_maledata2$Gene_ID

exp_maledata.2 = merge(traits, exp_maledata.1, by="row.names")

exp_maledata.2 = exp_maledata.2[
  with(exp_maledata.2, order(exp_maledata.2$Timepoint)),
]

hubSymbols = AnnotationDbi::mapIds(org.Mm.eg.db,
                                   keys = hubProbes_male,
                                   column = "SYMBOL",
                                   keytype = 'ENSEMBL') 

for(i in hubProbes_male) {
  y = exp_maledata.2[,i]
  sinreg_hub = sinreg(exp_maledata.2$Timepoint, y, plot=FALSE)
  
  y1 = exp_maledata.2[,i]
  y2 = sinreg_hub[[2]]
  
  ggplot(data=exp_maledata.2, aes(x=Timepoint, y=y1)) +
    geom_point() +
    geom_line(data=plotpoints, aes(x=Timepoint, y=y2, color="red")) +
    scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21)) +
    ggtitle(paste(names(hubProbes_male[which(hubProbes_male==i)]), "module hub gene: ", hubSymbols[which(hubProbes_male == i)])) +
    ylab("Expression") +
    xlab("Zeitgeber Time (ZT)") +
    theme_classic() +
    theme(legend.position="none")
  ggsave(paste("Sine_regression_plot_", ensembl[which(hubProbes_male == i)], "_ME", names(hubProbes_male[which(hubProbes_male==i)]), "_WT_clams.pdf", sep=""))

}

for(i in hubProbes_male) {
  y = exp_maledata.2[,i]
  sinreg_hub = sinreg(exp_maledata.2$Timepoint, y, plot=FALSE)
  
  y1 = exp_maledata.2[,i]
  y2 = sinreg_hub[[2]]
  
  ggplot(data=exp_maledata.2, aes(x=Timepoint, y=y1)) +
    geom_point() +
    geom_line(data=plotpoints, aes(x=Timepoint, y=y2, color="red")) +
    scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21)) +
    ggtitle(paste(names(hubProbes_male[which(hubProbes_male==i)]), "module hub gene: ", hubSymbols[which(hubProbes_male == i)])) +
    ylab("Expression") +
    xlab("Zeitgeber Time (ZT)") +
    theme_classic() +
    theme(legend.position="none")
  ggsave(paste("Sine_regression_plot_", hubSymbols[which(hubProbes_male == i)], "_ME", names(hubProbes_male[which(hubProbes_male==i)]), "_WT_clams.pdf", sep=""))
  
}
#### Heatmaps of R2 with FDRs for associations between Modules and Time using sinreg() ####

FDRs = as.matrix(results.ordered$FDR)
pVal = as.matrix(results.ordered$p.value)
rownames(pVal) = rownames(results.ordered)
colnames(pVal) = "pVal"

R2 = as.matrix(results.ordered$R2)
rownames(R2) = rownames(results.ordered)
colnames(R2) = "R2"

textMatrix = paste(ifelse((signif(FDRs, 1))<0.05, "*", ""), sep="")

tiff("Heatmap_WTfm_clams.tiff", res=400, height=7, width=2.5, units="in")
map1 = labeledHeatmap(Matrix = pVal,
                      xLabels = colnames(pVal),
                      yLabels = gsub("ME", "", rownames(pVal)),
                      ySymbols = rownames(pVal),
                      colorLabels = FALSE,
                      colors = brewer.pal(9, "RdBu"),
                      textMatrix = textMatrix,
                      setStdMargins = FALSE,
                      invertColors = FALSE,
                      cex.text = 1,
                      cex.lab.x = 0.75,
                      cex.lab.y = 0.65,
                      main=paste("Circadian Cycling of
Modules in Males and Females"))
dev.off()


save(results.ordered, FDRs, pVal, plotpoints, file=glue::glue("sinreg_association_with_gene_modules_WT.RData"))

#rain
rain_output <- rain(MEs_male, deltat=3, period=24, measure.sequence = c(12, 12, 12, 12, 12, 12, 10, 10), peak.border=c(0.3, 0.7), verbose=FALSE)
rain_output$FDR = p.adjust(rain_output$pVal, method="fdr")

mrain.ordered = rain_output[
  with(rain_output, order(rain_output$pVal, decreasing=FALSE)),
]
rownames(mrain.ordered) <- gsub(pattern = "ME", replacement = "", x = rownames(mrain.ordered), fixed = TRUE)
mrain.ordered$Module = rownames(mrain.ordered)

openxlsx::write.xlsx(mrain.ordered, file="rain_regresion_results_WT_clams.xlsx")

#rain heatmap
FDRs = as.matrix(mrain.ordered$FDR)
rownames(FDRs) = rownames(mrain.ordered)
colnames(FDRs) = "FDR"

pVal = as.matrix(mrain.ordered$pVal)
rownames(pVal) = rownames(mrain.ordered)
colnames(pVal) = "pVal"

textMatrix = paste(ifelse((signif(FDRs, 1))<0.05, "*", ""), sep="")

tiff("Heatmap_WTmf.rain.tiff", res=400, height=7, width=2.5, units="in")
map1 = labeledHeatmap(Matrix = pVal,
                      xLabels = colnames(pVal),
                      yLabels = gsub("ME", "", rownames(pVal)),
                      ySymbols = rownames(pVal),
                      colorLabels = FALSE,
                      colors = brewer.pal(9, "RdBu"),
                      textMatrix = textMatrix,
                      setStdMargins = FALSE,
                      invertColors = FALSE,
                      cex.text = 1,
                      cex.lab.x = 0.75,
                      cex.lab.y = 0.65,
                      main=paste("RAIN Rhythmic WT
Modules"))
dev.off()

#Female samples

MEnames = colnames(MEs_female2)
no.MEs = length(MEnames)

traits = pheno2$female$data

alldata = merge(traits, MEs_female2, by="row.names")

alldata = alldata[
  with(alldata, order(alldata$EntrainmentScore)),
]

head(alldata)
write.csv(alldata, "rw.ev.csv")
fitlist = as.list(1:no.MEs)
names(fitlist) <- MEnames

p.values = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(p.values) = MEnames
colnames(p.values) = "p.value"

R2 = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(R2) = MEnames
colnames(R2) = "R2"

period = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(period) = MEnames
colnames(period) = "period"

peak = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(peak) = MEnames
colnames(peak) = "peak"

amplitude = data.frame(matrix(ncol=1, nrow=no.MEs))
rownames(amplitude) = MEnames
colnames(amplitude) = "amplitude"

plotpoints = data.frame(matrix(ncol=no.MEs, nrow=nrow(alldata)))
rownames(plotpoints) = alldata$Row.names
colnames(plotpoints) = MEnames

for(i in MEnames){
  
  # print status
  print(paste("Running entity:", i, "which is", which(MEnames==i), "out of", no.MEs))
  
  #create temporary data matrix and model formula
  
  x = alldata$Timepoint
  y = alldata[,i]
  fitlist[[i]] <- sinreg(x, y, plot=FALSE)
  p.values[i,1] <- fitlist[[i]][[1]][6]
  R2[i,1] <- fitlist[[i]][[1]][5]
  period[i,1] <- fitlist[[i]][[1]][3]
  peak[i,1] <- fitlist[[i]][[1]][4]
  amplitude[i,1] <- fitlist[[i]][[1]][2]
  plotpoints[,i] <- fitlist[[i]][[2]]
}

results = cbind(p.values, R2, period, peak, amplitude)

results$FDR = p.adjust(results$p.value, method="fdr")
plotpoints$Timepoint = alldata$Timepoint

for(i in MEnames) {
  y = alldata[,i]
  y1 = plotpoints[,i]
  ggplot(data=alldata, aes(x=Timepoint, y=y)) +
    geom_point() +
    geom_line(data=plotpoints, aes(x=Timepoint, y=y1, color="red")) +
    scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21)) +
    ggtitle(paste(i)) +
    ylab("Module EigenValue") +
    xlab("Zeitgeber Time (ZT)") +
    theme_classic() +
    theme(legend.position="none")
  ggsave(paste("Sine_regression_plot_", i, "_females.pdf"))
}

# order results by best fit (highest R2 value) and save results #
results.ordered = results[
  with(results, order(results$R2, decreasing=TRUE)),
]

results.ordered$Module = rownames(results.ordered)

openxlsx::write.xlsx(results.ordered, file="Sine_regresion_results_females.xlsx")

# Checking if hub genes within modules cycle in a similar manner to EigenValues #

exp_femdata.1 = as.data.frame(t(exp_femdata2[,-1]))
colnames(exp_femdata.1) = exp_femdata2$Gene_ID

exp_femdata.2 = merge(traits, exp_femdata.1, by="row.names")

exp_femdata.2 = exp_femdata.2[
  with(exp_femdata.2, order(exp_femdata.2$GenotypeScores)),
]  
#exp_femdata.2 = exp_femdata.2[
#  with(exp_femdata.2, order(exp_femdata.2$Timepoint)),
#]

hubSymbols = AnnotationDbi::mapIds(org.Mm.eg.db,
                                   keys = hubProbes_female,
                                   column = "SYMBOL",
                                   keytype = 'ENSEMBL') 

for(i in hubProbes_female) {
  y = exp_femdata.2[,i]
  sinreg_hub = sinreg(exp_femdata.2$Timepoint, y, plot=FALSE)
  
  y1 = exp_femdata.2[,i]
  y2 = sinreg_hub[[2]]
  
  ggplot(data=exp_femdata.2, aes(x=Timepoint, y=y1)) +
    geom_point() +
    geom_line(data=plotpoints, aes(x=Timepoint, y=y2, color="red")) +
    scale_x_continuous(breaks=c(0, 3, 6, 9, 12, 15, 18, 21)) +
    ggtitle(paste(names(hubProbes_female[which(hubProbes_female==i)]), "module hub gene: ", hubSymbols[which(hubProbes_female == i)])) +
    ylab("Expression") +
    xlab("Zeitgeber Time (ZT)") +
    theme_classic() +
    theme(legend.position="none")
  ggsave(paste("Sine_regression_plot_", hubSymbols[which(hubProbes_female == i)], "_ME", names(hubProbes_female[which(hubProbes_female==i)]), "_females.pdf", sep=""))
  
}

#### Heatmaps of R2 with FDRs for associations between Modules and Time using sinreg() ####

FDRs = as.matrix(results.ordered$FDR)
pVal = as.matrix(results.ordered$p.value)
rownames(pVal) = rownames(results.ordered)
colnames(pVal) = "pVal"

R2 = as.matrix(results.ordered$R2)
rownames(R2) = rownames(results.ordered)
colnames(R2) = "R2"

textMatrix = paste(ifelse((signif(FDRs, 1))<0.05, "*", ""), sep="")

tiff("Heatmap_females.tiff", res=400, height=7, width=2.5, units="in")
map1 = labeledHeatmap(Matrix = pVal,
                      xLabels = colnames(pVal),
                      yLabels = gsub("ME", "", rownames(pVal)),
                      ySymbols = rownames(pVal),
                      colorLabels = FALSE,
                      colors = brewer.pal(9, "RdBu"),
                      textMatrix = textMatrix,
                      setStdMargins = FALSE,
                      invertColors = FALSE,
                      cex.text = 1,
                      cex.lab.x = 0.75,
                      cex.lab.y = 0.65,
                      main=paste("Circadian Cycling of
Modules in females"))
dev.off()


save(results.ordered, FDRs, pVal, plotpoints, file=glue::glue("sinreg_association_with_gene_modules_females.RData"))

#rain
rain_output <- rain(MEs_female, deltat=3, period=24, measure.sequence = c(6, 6, 6, 6, 6, 6, 5, 5), peak.border=c(0.3, 0.7), verbose=FALSE)
rain_output$FDR = p.adjust(rain_output$pVal, method="fdr")

frain.ordered = rain_output[
  with(rain_output, order(rain_output$FDR, decreasing=FALSE)),
]
rownames(frain.ordered) <- gsub(pattern = "ME", replacement = "", x = rownames(frain.ordered), fixed = TRUE)
frain.ordered$Module = rownames(frain.ordered)

openxlsx::write.xlsx(frain.ordered, file="rain_regresion_results_females.xlsx")

#rain heatmap
FDRs = as.matrix(frain.ordered$FDR)
rownames(FDRs) = rownames(frain.ordered)
colnames(FDRs) = "FDR"

pVal = as.matrix(frain.ordered$pVal)
rownames(pVal) = rownames(frain.ordered)
colnames(pVal) = "pVal"

textMatrix = paste(ifelse((signif(FDRs, 1))<0.05, "*", ""), sep="")

tiff("Heatmap_females.rain.tiff", res=400, height=7, width=2.5, units="in")
map1 = labeledHeatmap(Matrix = pVal,
                      xLabels = colnames(pVal),
                      yLabels = gsub("ME", "", rownames(pVal)),
                      ySymbols = rownames(pVal),
                      colorLabels = FALSE,
                      colors = brewer.pal(9, "RdBu"),
                      textMatrix = textMatrix,
                      setStdMargins = FALSE,
                      invertColors = FALSE,
                      cex.text = 1,
                      cex.lab.x = 0.75,
                      cex.lab.y = 0.65,
                      main=paste("RAIN Rhythmic
Modules in females"))
dev.off()

#name extractor and gene symbol generator
male_modulemem <- read.delim("MALES Probe Module Membership.txt", sep = "\t",
                             header = TRUE, stringsAsFactors = FALSE)
male_ens.mods <- setDT(male_modulemem[,c(11:12)])
#test <- DT[order(Module)]
#test <- DT[Module == "black"]

female_modulemem <- read.delim("FEMALES Probe Module Membership.txt", sep = "\t",
                               header = TRUE, stringsAsFactors = FALSE)
female_ens.mods <- setDT(female_modulemem[,c(11:12)])

msamples.vector <- mrain.ordered$Module
modules_interest = msamples.vector
for (i in modules_interest) {
  
  gene_id <- colnames(exp$maledata$data)[consensusMods$colors == i]
  ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  genes <- getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", 
                 values = gene_id, mart = ensembl) %>% unlist %>% as.character %>% unique %>% sort
  write.csv(genes, (paste(i,"_male_geneslist.csv")), sep = "\t", quote = FALSE)
}

fsamples.vector <- frain.ordered$Module
modules_interest = fsamples.vector
for (i in modules_interest) {
  
  gene_id <- colnames(exp$femdata$data)[fMods$colors == i]
  ensembl <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  genes <- getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", 
                 values = gene_id, mart = ensembl) %>% unlist %>% as.character %>% unique %>% sort
  write.csv(genes, (paste(i,"_female_geneslist.csv")), sep = "\t", quote = FALSE)
}

#print full list with genes#
gene_id_male <- colnames(exp$maledata$data)[consensusMods$colors == i]
#ensembl_male <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
genes_interest_male <- getBM(attributes = "external_gene_name", filters = "ensembl_gene_id", 
                             values = probes_interest_male, mart = ensembl_male) %>% unlist %>% as.character %>% unique %>% sort
write.csv(genes_interest_male, file = paste("module_interest_", modules_interest[which(modules_interest == i)], "_males.csv", sep=""))
}
###EnrichR Pathway Analysis###

modules_interest = c("green", "grey", "black", "red", "blue", "brown", "yellow", "pink", "magenta", "turquoise")
lapply(modules_interest, function(module) {
  data = read.csv(glue::glue("{module} _male_geneslist.csv")) 
  
  data %>%
    dplyr::select(x) %>%
    purrr::flatten() %>%
    enrichR::enrichr(c("GO_Biological_Process_2018",
                       "GO_Cellular_Component_2018",
                       "GO_Molecular_Function_2018",
                       "KEGG_2019_Mouse",
                       "Panther_2016",
                       "Reactome_2016",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>% 
    purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = "")) %T>%
    #purrr::map(~ dplyr::filter(., Adjusted.P.value < 0.05)) %>% 
    #purrr::map(~ dplyr::filter(., stringr::str_detect(Genes, ";"))) %>% 
    openxlsx::write.xlsx(file = glue::glue("Module_{module}_males_enrichr.xlsx")) %>%
    DMRichR::slimGO(tool = "enrichR",
                    annoDb = "org.Mm.eg.db",
                    plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("Module_{module}_males_rrvgo_enrichr.xlsx")) %>%
    DMRichR::GOplot() %>%
    ggplot2::ggsave(glue::glue("Module_{module}_males_enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10) 
  
})

## testing
male_modulemem <- read.delim("MALES Gene Module Membership.txt", sep = "\t",
                                header = TRUE, stringsAsFactors = FALSE)
DT <- setDT(male_modulemem[,c(21:22)])
test <- DT[order(Module)]
test <- DT[Module == "black"]

exp <- list(femdata = list(data = as.data.frame(t(exp_femdata[,-c(1)]))),
            maledata = list(data = as.data.frame(t(exp_maledata[,-c(1)]))))
names(exp$femdata$data) = exp_femdata$Gene_ID
rownames(exp$femdata$data) = names(exp_femdata)[-c(1)]
names(exp$maledata$data) = exp_maledata$Gene_ID
rownames(exp$maledata$data) = names(exp_maledata)[-c(1)]
checkSets(exp)

test_list <- list(modules = list(data = as.data.frame(t(male_modulemem[,c(21:22)]))))
names(test_list$modules$data) = male_modulemem$Module
rownames(test_list$modules$data) = names(DT)
names(exp$maledata$data) = exp_maledata$Gene_ID
rownames(exp$maledata$data) = names(exp_maledata)[-c(1)]
checkSets(test_list)

seto
test2 <- male_modulemem[with(male_modulemem, Module > pink), ]
test3 <- DT[, .N, by = Module]
test4 <- DT[Module == "pink", .N, by = Module]
vector_test <- c(mrain.ordered$Module)
names(vector_test)[names(vector_test) == "c(mrain.ordered$Module)"] <- "Module"
testset_modules <- data[match(vector_test, male_modulemem$Module)]
testset_modules <- order_by(vector_test, male_modulemem$Module)
mMods2 <- moduleMembership2$Module
MEs_male <- moduleEigengenes(t(exp_maledata[,-c(1)]), colors = mMods2)$eigengenes
rownames(MEs_male) <- rownames(t(exp_maledata[,-c(1)]))
mMEs <- list(male = list(data = MEs_male))
mMEs <- orderMEs(mMEs)
order_by()
female_modulemem <- read.delim("FEMALES Gene Module Membership.txt", sep = "\t",
                                header = TRUE, stringsAsFactors = FALSE)
fMods2 <- moduleMembership2$Module
MEs_female <- moduleEigengenes(t(exp_femdata[,-c(1)]), colors = fMods2)$eigengenes
rownames(MEs_female) <- rownames(t(exp_femdata[,-c(1)]))
fMEs <- list(male = list(data = MEs_female))
fMEs <- orderMEs(fMEs)

rm(moduleMembership2)


