
cf_trait <- read.csv("cf_traits.csv")
cf_trait <- read.csv("cf_traits_test.csv")
cm_trait <- read.csv("cm_traits.csv")
rf_trait <- read.csv("rf_traits.csv")
rm_trait <- read.csv("rm_traits.csv")

traits = as.data.frame(cf_trait)
traits = as.data.frame(rm_trait)
rownames(traits) <- traits$ID
traits = as.data.frame(pheno$cm$data)
traits = as.data.frame(pheno$wtf$data)
traits = as.data.frame(pheno$wtm$data)
rownames(traits) <- traits$SampleID
#row.names(traits) <- traits$sampleID
row.names(traits)
row.names(MEnames)
row.names(consensusMEs$rmdata$data)
alldata = merge(traits, MEs_cm, by="row.names")
alldata = merge(traits, MEs_wt, by="row.names")
alldata_cf = merge(traits, consensusMEs$cfdata$data, by="row.names")
alldata_rm = merge(traits, consensusMEs$rmdata$data, by="row.names")
alldata = merge(traits, consensusMEs$wtmdata$data, by="row.names")

alldata = alldata[
  with(alldata, order(alldata$Timepoint)),
]

ggplot(alldata_cf, aes(x = alldata_cf$WTvDEL, y = alldata_cf$MEsaddlebrown, fill = Genotype)) +
  geom_boxplot() +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text = element_text(face = "bold", color = "black")) +
  ylab("") +
  xlab("Genotypes") +
  ggtitle("Module")

library(ggpubr)
library(rstatix)
ggplot(alldata_cf, aes(x = alldata_cf$LD.Avg.Light.XTOT, y = alldata_rf$MEsaddlebrown, fill = Genotype)) +
  #coord_flip() +
  geom_boxplot() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text = element_text(face = "bold", color = "black")) +
  ylab("ModuleEigenvalue") +
  xlab("Genotypes RF") +
  ggtitle("Module Turquoise")
ggboxplot(alldata_rf, x = "Genotype", y = "MEturquoise", fill = "Genotype")
ggplot(alldata_cf, aes(x = alldata_cf$LD.Avg.Light.XTOT, y = alldata_cf$MEsaddlebrown, color = Genotype)) +
  geom_point(size = 10) +
#  stat_ellipse(type = "norm", level = 0.5) +
#  scale_color_gradient(low = "red", high = "blue") +
#  geom_smooth(method = "loess") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text = element_text(face = "bold", color = "black")) +
  ylab("ModuleEigenvalue") +
  xlab("LD Avg Light XTOT") +
  ggtitle("Module Saddlebrow")


g4 + stat_compare_means(method = "anova")
my_comparisons <- list(c("WTWT11", "WTWT12"), c("HETWT11", "HETWT12"), c("WTTG11", "WTTG12"), c("HETTG11", "HETTG12"), c("WTWT11", "HETWT11"), c("WTWT12", "HETWT12"), c("WTWT11", "WTTG11"), c("WTWT12", "WTTG12"), c("WTWT11", "HETTG11"), c("WTWT12", "HETTG12"))
my_comparisons <- list(c("WTWT11", "HETTG11"), c("WTWT12", "HETTG12"), c("WTWT11", "HETWT11"), c("WTWT12", "HETWT12"), c("WTWT11", "WTTG11"), c("WTWT12", "WTTG12"))
my_comparisons <- list(c("WTWT11", "HETWT11"), c("WTWT12", "HETWT12"))
g3 <- ggboxplot(alldata_cf, x = "WTvDEL", y = "MEsaddlebrown", fill = "WTvDEL")
g4 + stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = TRUE)
test <- alldata_cf
test$MEsaddlebrown <- as.factor(test$MEsaddlebrown)
test$Genotype <- as.(test$Genotype)
stat.test <- aov(Entrainment ~ MEsaddlebrown, data = test) %>% tukey_hsd()
stat.test <- tutest %>% add_xy_position(x = "Genotype", fun = "mean_sd", dodge = 0.8)

ggboxplot(alldata_cf, x = "Genotype", y = "MEsaddlebrown") 
  stat_pvalue_manual(
    stat.test, label = "p.adj",
    y.position = c(.32, .34, .36, .38, .40, .42)
  )
p + stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format", paired = TRUE)
p
alldata_cf <- alldata_cf[order(alldata_cf$Genotype, decreasing = TRUE),]
bxp <- ggboxplot(alldata_cf, x = "Genotype", y = "MEsaddlebrown", fill = "Genotype")
bxp

alldata_cf = merge(traits, consensusMEs$cfdata$data, by="row.names")

alldata_cf$GenotypeScores <- as.factor(alldata_cf$Genotype)
alldata_cf$GenotypeScores <- as.numeric(alldata_cf$GenotypeScores)
alldata_cf$Entrainment <- as.factor(alldata_cf$Entrainment)

tutest <- aov(Entrainment ~ MEsaddlebrown*Genotype, data = alldata_cf) %>% tukey_hsd()
alldata_cf$MEsaddlebrown <- as.factor(alldata_cf$MEsaddlebrown)
tutest <- aov(MEsaddlebrown ~ Entrainment*Genotype, data = alldata_cf) %>% tukey_hsd()

bxp + stat_pvalue_manual(tutest, label = "p.adj", y.position = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14))
tutest <- tutest %>% add_xy_position(x = "Genotype", fun = "mean_sd", dodge = 0.8)
alldata_cf$MEsaddlebrown <- as.factor(alldata_cf$MEsaddlebrown)
alldata_cf$MEturquoise <- as.factor(alldata_cf$MEturquoise)

test4 <- alldata_cf %>% group_by(Genotype, Entrainment) #%>% summarize(mean_saddlebrown = mean(MEsaddlebrown))
head(test4)
test5 <- test4 %>% summarise(sb = mean(MEsaddlebrown), tq = mean(MEturquoise))
model <- aov(MEsaddlebrown ~ Genotype*Entrainment, data = test4)
summary(model)
