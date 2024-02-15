library(WGCNA)
BiocManager::install("gghalves")
BiocManager::install("devtools")

devtools::install_github('smin95/smplot2', force = TRUE)

library(gghalves)
library(devtools)
library(smplot2)
load("MEMAs.RData")
getwd()

head(MEMAs)
cor_ENT <- data.frame(MEMAs[[1]]$ID, MEMAs[[1]]$corPearson.female, MEMAs[[1]]$p.equalWeights)
cor_ENT$ID = cor_ENT$MEMAs..1...ID
cor_Gen <- data.frame(MEMAs[[2]]$ID, MEMAs[[2]]$cor.female, MEMAs[[2]]$p.equalWeights)
cor_Gen$ID = cor_Gen$MEMAs..2...ID
cor_Sex <- data.frame(MEMAs[[3]]$ID, MEMAs[[3]]$corPearson.female, MEMAs[[3]]$p.equalWeights)
cor_Sex$ID = cor_Sex$MEMAs..3...ID
cor_Merged <- merge(cor_ENT, cor_Gen, by = "ID")
cor_Merged <- merge(cor_Merged, cor_Sex, by = "ID")
head(cor_Merged)

cor_Merged$MEMAs..1...ID <- NULL
cor_Merged$MEMAs..2...ID <- NULL
cor_Merged$MEMAs..3...ID <- NULL
cor_maxtrix <- as.matrix(cor_Merged)
head(cor_maxtrix)

star <- apply(pvalues, 2, function(x){sapply(x, function(y){ifelse(y < 0.05, "*", "")})})
pdf("Consensus Modules Meta Covariate Correlation Plot CLAMS-RW3.pdf", width = 11, height = 15)
sizeGrWindow(width = 11, height = 15)
par(mar = c(9, 8, 1, 2))
labeledHeatmap(Matrix = cor_maxtrix, xLabels = colnames(cor_Merged), yLabels = rownames(zscores), 
               ySymbols = gsub("ME", "", rownames(zscores)), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = star, setStdMargins = FALSE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-5, 5), main = "", cex.lab.y = 1)
dev.off()

pheno2 <- read.csv("phospho.ev.csv")
pheno2 <- pheno2[,-c(1)]
# Set up variables to contain the module-trait correlations
moduleTraitCor = list();
moduleTraitPvalue = list();
# Calculate the correlations
for (set in 1:nSets)
{
  moduleTraitCor = cor(consensusMEs$female$data, pheno2$female$data, method =  c("pearson", "kendall", "spearman"), use = "p");
  moduleTraitPvalue = corPvalueFisher(moduleTraitCor, 124);
}

moduleTraitCorP = cor(consensusMEs$female$data, pheno2$female$data, method = "pearson", use = "p")
moduleTraitPvalueP2 = corPvalueFisher(moduleTraitCorP, 30, twoSided = TRUE)

moduleTraitCorP = cor(MEs_female2, cov_female, method = "pearson", use = "p")
moduleTraitPvalueP2 = corPvalueFisher(moduleTraitCorP, 31, twoSided = TRUE)
#moduleTraitQvalueP = as.matrix(p.adjust(moduleTraitPvalueP2, method = p.adjust.methods, n = length(moduleTraitPvalueP2)))

moduleTraitCorS = cor(consensusMEs$female$data, pheno2$female$data, method = "spearman", use = "p")
moduleTraitPvalueS2 = corPvalueFisher(moduleTraitCorS, 30, twoSided = TRUE)
#moduleTraitQvalueS = p.adjust(p, method = p.adjust.methods, n = length(p))

moduleTraitCorK = cor(consensusMEs$female$data, pheno2$female$data, y = NULL, method = "kendall", use = "p")
moduleTraitPvalueK2 = corPvalueFisher(moduleTraitCorK, 32, twoSided = TRUE)
#moduleTraitQvalueK = p.adjust(p, method = p.adjust.methods, n = length(p))

library(ggpubr)
test <- read.csv("clams_cf.me.traits.csv")
test <- read.csv("clams_cf.me.traits.csv")
test <- read.csv("rw_rm.me.traits.csv")
test <- read.csv("test2.csv")
ggscatter(test, x = "AvgKcalIntake", y = "MEblack", color = "Genotype", shape = "Genotype",
          pallete = c("#00AFBB", "#E7B800", "#FC4E07", "#CC79A7"),
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Entrainment", ylab = "module eigenvalue")

ggscatter(test, x = "AvgKcalIntake", y = "MEblack", color = "genotype",
          pallete = c("#00AFBB", "#E7B800", "#FC4E07", "#CC79A7"))

ggscatter(df, x = "AvgKcalIntake", y = "MEblack",
          color = "Genotype", shape = "Genotype",
          pallete = c("#00AFBB", "#E7B800", "#FC4E07", "#CC79A7"),
          ellipse = TRUE, mean.point = TRUE,
          star.plot = TRUE)

cor.test(test$PeriodFirst4, test$MEdarkred, method = "spearman")
cormat <- round(cor(test, use = 'pairwise.complete.obs', method = "spearman"),2)
head(cormat)
library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)
library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()
tmp <- cor(test, test, method = "spearman")

ggplot(data = test, mapping = aes(x = AvgKcalIntake, y = MEblack)) +
  geom_point(shape = 21, fill = '#0f993d', color = 'red', size = 3) +
  sm_statCorr(color = '#0f993d', corr_method = 'spearman',
              linetype = 'dashed')

ggplot(data = test, mapping = aes(x = AvgKcalIntake, y = MEblack, color = genotype)) + 
  geom_point(size = 2)
unique_genotypes <- unique(test$genotype)
number_of_genotypes <- length(unique_genotypes)
number_of_genotypes

ggplot(data = test, mapping = aes(x = AvgKcalIntake, y = MEblack, color = genotype)) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("red", "green", "yellow", "blue")) + 
  sm_corr_theme()

ggplot(data = test, mapping = aes(x = AvgKcalIntake, y = MEblack, color = factor(genotype))) +
  geom_point(size = 2) + 
  scale_color_manual(values = c("red", "green", "yellow", "blue"))

test <- read.csv("rw_ev.traits.csv")
test <- read.csv("cf_ev.traits.csv")
test <- read.csv("cm_ev.traits.csv")

rw_ev_traits <- read.csv("rwm_eigenvalue_organized.csv")
cf_ev_traits <- read.csv("cf_eigenvalue_organized.csv")
cm_ev_traits <- read.csv("cm_eigenvalue_organized.csv")
names(MEs_cf)
row.names(MEs_cf)

c
# Use this to extract individual correlations/genotype and entrainment group
ggplot(data = rw_ev_traits, mapping = aes(x = PeriodBasal, y = MEpink, color = factor(Genotype))) +
  ggtitle("RM ME Black vs PeriodBasal Spearman Correlation") +
  geom_point() +
  #scale_color_manual(values = c("red", "green", "yellow", "blue")) +
  #labs(color = factor("genotype")) +
  #theme(legend.position = "bottom") +
  stat_ellipse(type = "norm", level = 0.5) +
  #sm_statCorr(color = 'black', corr_method = 'spearman', linetype = 'dashed') +
  sm_corr_theme(borders = FALSE, legends = TRUE)

ggplot(data = rw_ev_traits, mapping = aes(x = c("PeriodBasal", "PeriodFirst4", "PeriodLast4", "AmplitudeBasal", "AmplitudeFirst4", "AmplitudeLast4"), y = names(MEs_rw), color = factor(Genotype))) +
  ggtitle("RM ME Black vs PeriodBasal Spearman Correlation") +
  geom_point() +
  #scale_color_manual(values = c("red", "green", "yellow", "blue")) +
  #labs(color = factor("genotype")) +
  #theme(legend.position = "bottom") +
  stat_ellipse(type = "norm", level = 0.5) +
  #sm_statCorr(color = 'black', corr_method = 'spearman', linetype = 'dashed') +
  sm_corr_theme(borders = FALSE, legends = TRUE)

# Use this to generate corrlation plot
ggplot(data = rw_ev_traits, mapping = aes(x = PeriodBasal, y = MEpink)) +
  ggtitle("RM ME Pink vs PeriodBasal Spearman Correlation") +
  geom_point() +
  sm_statCorr(color = 'black', corr_method = 'spearman',
              linetype = 'dashed')

###Did not work###
plot_cor <- function(u,v){
  ggplot(data = rw_ev_traits, mapping = aes(x = u, y = v)) +
    ggtitle(glue::glue("RM ME {v} vs {u} Spearman Correlation")) +
    geom_point() +
    sm_statCorr(color = 'black', corr_method = 'spearman',
                linetype = 'dashed')
  ggsave(glue::glue("RW ME {v} vs {u}.pdf"))
}
mapply(plot_cor,test_list2,test_list1)

###Forloop version###
###Running wheel test list###
rw_traits <- c("PeriodBasal", "PeriodFirst4", "PeriodLast4", "AmplitudeBasal", "AmplitudeFirst4", "AmplitudeLast4")

###CLAMS test list###
clams_traits <- c("AvgRER",	"AvgLightRER",	"AvgDarkRER",	"AvgHeat",	"AvgLightHeat",	"AvgDarkHeat",	"AvgKcalIntake", "AvgLightKcalIntake",	"AvgDarkKcalIntake")
###Set-up data###
test_list1 <- names(MEs_cm)
test_list2 <- clams_traits
data = cm_ev_traits
p = as.data.frame(moduleTraitPvaluecm)

###Creates correlation plots###

for (i in test_list2) {
  for (j in test_list1) {

x1 = data[,i]
y1 = data[,j]

if (p[j,i] < 0.05) {
  ggplot(data, mapping = aes(x = x1, y = y1)) +
    ggtitle(paste("RW", j, "vs", i, "Spearman Correlation")) +
    geom_point() +
    sm_statCorr(color = 'black', corr_method = 'spearman',
                linetype = 'dashed')
  ggsave(paste0("RW_", j, "_vs_", i, ".pdf"))
}
else {
  print("No significant data")
}}}

###Creates ellipses###
for (i in test_list2) {
  for (j in test_list1) {
    
    x1 = data[,i]
    y1 = data[,j]
    
    if (p[j,i] < 0.05) {
      ggplot(data, mapping = aes(x = x1, y = y1, color = factor(Genotype))) +
        ggtitle(paste("CM", j, "vs", i, "Spearman Correlation")) +
        geom_point() +
        stat_ellipse(type = "norm", level = 0.5) 
      ggsave(paste0("CM_", j, "_vs_", i, "ellipse.pdf"))
    }
    else {
      print("No significant data")
    }}}

