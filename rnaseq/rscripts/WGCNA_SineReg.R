#### Using sine regression to assess associations between WGCNA gene modules and time ####

library(WGCNA)
library(ShellChron)
library(ggplot2)
library(openxlsx)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(RColorBrewer)

setwd("/Users/karineier/Documents/Circadian-Analysis/Circadian-Analysis")

load("WGNCAoptimization6.8.21.RData")


MEnames = colnames(MEs_female)
no.MEs = length(MEnames)

traits = pheno$female$data

alldata = merge(traits, MEs_female, by="row.names")

alldata = alldata[
  with(alldata, order(alldata$Timepoint)),
]

head(alldata)

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
  
  #assign fit to list by name
  
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
  ggsave(paste("Sine_regression_plot_", i, "_Females.pdf"))
}

# order results by best fit (highest R2 value) and save results #
results.ordered = results[
  with(results, order(results$R2, decreasing=TRUE)),
]

results.ordered$Module = rownames(results.ordered)

openxlsx::write.xlsx(results.ordered, file="Sine_regresion_results_females.xlsx")

# Checking if hub genes within modules cycle in a similar manner to EigenValues #

exp_femdata.1 = as.data.frame(t(exp_femdata[,-1]))
colnames(exp_femdata.1) = exp_femdata$Gene_ID

exp_femdata.2 = merge(traits, exp_femdata.1, by="row.names")

exp_femdata.2 = exp_femdata.2[
  with(exp_femdata.2, order(exp_femdata.2$Timepoint)),
]

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
  ggsave(paste("Sine_regression_plot_", hubSymbols[which(hubProbes_female == i)], "_ME", names(hubProbes_female[which(hubProbes_female==i)]), "_Females.pdf", sep=""))
  
}

#### Heatmaps of R2 with FDRs for associations between Modules and Time using sinreg() ####

FDRs = as.matrix(results.ordered$FDR)
R2 = as.matrix(results.ordered$R2)
rownames(R2) = rownames(results.ordered)
colnames(R2) = "R2"

textMatrix = paste(ifelse((signif(FDRs, 1))<0.05, "*", ""), sep="")

tiff("Heatmap_Females.tiff", res=400, height=7, width=2.5, units="in")
map1 = labeledHeatmap(Matrix = R2,
                      xLabels = colnames(R2),
                      yLabels = gsub("ME", "", rownames(R2)),
                      ySymbols = rownames(R2),
                      colorLabels = FALSE,
                      colors = brewer.pal(9, "OrRd"),
                      textMatrix = textMatrix,
                      setStdMargins = FALSE,
                      cex.text = 1,
                      cex.lab.x = 0.75,
                      cex.lab.y = 0.65,
                      main=paste("Circadian Cycling of
Modules in Females"))
dev.off()


save(results.ordered, FDRs, R2, plotpoints, file=glue::glue("sinreg_association_with_gene_modules_females.RData"))


