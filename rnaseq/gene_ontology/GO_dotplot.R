#Author : Osman Sharifi

####################
## load libraries ##
####################
library(dplyr)
library(ggplot2)
library(readxl)
BiocManager::install("rio")
library(rio)
library(viridis)

###############
## load data ##
###############
module <- import_list("/Users/aron/Desktop/LaSalle_Lab/Analysis/clamsrw/rnaseq/separate_normalization/GOTermsUpdatedCorrected/CFGOTerms/Module_black_CF8_enrichr.xlsx", rbind = TRUE)
colnames(module) <- c("Term", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value", "Old.Adjusted.Pl.value", "Odds.Ratio", "Combined.Score", "Genes", "Ontology")
module <- module %>%
  select(Term, Overlap, Adjusted.P.value, Genes, Odds.Ratio, Ontology) %>%
  filter(Adjusted.P.value <= 0.05)
module <- module %>% 
  mutate(Ontology = ifelse(as.character(Ontology) == 1, "BP" , as.character(Ontology)))
module <- module %>% 
  mutate(Ontology = ifelse(as.character(Ontology) == 2, "CC" , as.character(Ontology)))
module <- module %>% 
  mutate(Ontology = ifelse(as.character(Ontology) == 3, "MF" , as.character(Ontology)))
df = module[module$Ontology == "BP" | module$Ontology == "MF" | module$Ontology == "CC", ]

#make unique order based on Odds.Ratio
df.order <- unique(as.character(df$Term)[order(df$Odds.Ratio, decreasing = FALSE)])

#reassign journal as factor with new levels
df$Term <- factor(df$Term, levels = df.order)

########################
## plot filtered data ##
########################
ggplot(df,
       aes(x = Ontology, y = Term , size = Odds.Ratio, fill = Adjusted.P.value)) +
  geom_point(shape = 21) +
  scale_fill_viridis() + 
  xlab('') + ylab('') +
  labs(
    title = 'Top enrichR Terms',
    subtitle = 'Green Module '
  ) 

