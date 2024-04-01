### Pathway Analaysis for Select Modules ###

setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/WGCNAoptimization/CircadianModules/Males/GOterms")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/WGCNA_Final/circadianmodules/GoTerms")
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/clamsrw/rnaseq/separate_normalization/GOTermsUpdatedCorrected")
packages <- c("edgeR", "tidyverse", "magrittr", "RColorBrewer", "org.Mm.eg.db", "AnnotationDbi",
              "enrichR", "openxlsx", "gt", "glue", "DMRichR")
BiocManager::install("enrichR")
BiocManager::install("DMRichR")
BiocManager::install("clusterProfiler")
BiocManager::install("glue")
install.packages("DMRichR")
library(enrichR)
library(DMRichR)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(magrittr)
library(glue)
library(plyr)
library(data.table)

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = TRUE)
BiocManager::install("ben-laufer/DMRichR")
if(R.Version()$major < 4)
  install.packages("dmrseq", repos = "https://bioconductor.org/packages/3.12/bioc")

enrichR:::.onAttach() # Needed or else "EnrichR website not responding"
stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))

### Loading Data ###

modules_interest = c("pink", "purple", "red", "salmon", "darkturquoise")
modules_interest = c("yellow", "pink", "black")
lapply(modules_interest$Module, function(module) {
  data = read.csv(glue::glue("{module}_universal_modules.csv")) 
  
  data %>%
    dplyr::select(x) %>%
    purrr::flatten() %>%
    enrichR::enrichr(c("GO_Biological_Process_2023",
                       "GO_Cellular_Component_2023",
                       "GO_Molecular_Function_2023",
                       "KEGG_2019_Mouse",
                       "dbGaP",
                       "Panther_2016",
                       "Reactome_2022",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>% 
    purrr::set_names(names(.) %>% stringr::str_trunc(31, ellipsis = "")) %T>%
    #purrr::map(~ dplyr::filter(., Adjusted.P.value < 0.05)) %>% 
    #purrr::map(~ dplyr::filter(., stringr::str_detect(Genes, ";"))) %>% 
    openxlsx::write.xlsx(file = glue::glue("Module_{module}_Universal_enrichr.xlsx")) 
})

###generate GO Term list###
test_modules <- lapply(modules_interest$Module, function(module) {
  data = read.csv(glue::glue("{module}_universal_modules.csv")) 

  data %>%
    dplyr::select(x) %>%
    purrr::flatten() %>%
    enrichR::enrichr(c("GO_Biological_Process_2023",
                       "GO_Cellular_Component_2023",
                       "GO_Molecular_Function_2023",
                       "KEGG_2019_Mouse",
                       "dbGaP",
                       "Panther_2016",
                       "Reactome_2022",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>% 
    purrr::set_names(names(.) %>% stringr::str_trunc(40, ellipsis = "")) %T>%
    purrr::map(~ dplyr::filter(., Adjusted.P.value < 0.05))
    #new_name <- paste0("test_RW",module) %>%
    #assign(new_name, test_modules)
    #purrr::map(~ dplyr::filter(., stringr::str_detect(Genes, ";"))) %>% 
    #openxlsx::write.xlsx(file = glue::glue("Module_{module}_RM12_enrichr.xlsx")) 
})

###names list items by module
names(test_modules) <- modules_interest$Module

##extract KEGG##
Kegg <- lapply(modules_interest$Module, function(module) {
  test_modules[[module]]$KEGG_2019_Mouse  
  })
names(Kegg) <- modules_interest$Module

##bind KEGG##
###create sample column and name by module color###
##for loop worked!###
for(i in modules_interest$Module) {
  Kegg[[i]]$Sample <- paste0(i)  
}
Kegg[["darkslateblue"]]$Sample <- paste0("darkslateblue")  
###Binds all list items###
Kegg_bound <- rbindlist(Kegg, fill = TRUE)
Kegg_filtered <- filter(Kegg_bound, Adjusted.P.value <= 0.05)

##extract KEGG##
trial2 <- lapply(modules_interestrm, function(module) {
  test_modules_rw[[module]]$KEGG_2019_Mouse  
})
names(trial2) <- modules_interestrm

##bind KEGG##
###creates sample column and name by module color###
##for loop worked!###
for(i in modules_interestrm) {
  trial2[[i]]$Sample <- paste0("RW_", i)  
}

###Binds all list items###
test11 <- rbindlist(trial2)
test11 <- filter(test11, Adjusted.P.value <= 0.05)
Kegg_rm_bound <- test11
All_exp <- rbind(Kegg_rm_bound, Kegg_cf_bound, Kegg_cm_bound)
DD_exp <- rbind(Kegg_cf_bound, Kegg_cm_bound)
Male_exps <- rbind(Kegg_rm_bound, Kegg_cm_bound)

go_analysis <- read.delim(read.csv("grey_moduleCM16.csv"))
simMatrix <- calculateSimMatrix(test_modules_cf$turquoise$GO_Biological_Process_2023,
                                orgdb = "org.Mm.eg.db",
                                ont = "BP",
                                method = "Rel")

#make unique order based on Odds.Ratio
Kegg_cf_bound.order <- unique(as.character(Kegg_cf_bound$Term)[order(Kegg_cf_bound$Adjusted.P.value, decreasing = TRUE)])
kegg_rhythmic.order <- unique(as.character(kegg_rhythmic$Term)[order(kegg_rhythmic$Adjusted.P.value, decreasing = TRUE)])
All_exp.order <- unique(as.character(All_exp$Term)[order(All_exp$Adjusted.P.value, decreasing = TRUE)])
Male_exps.order <- unique(as.character(Male_exps$Term)[order(Male_exps$Adjusted.P.value, decreasing = TRUE)])
#reassign journal as factor with new levels
Kegg_rm_bound2 <- Kegg_rm_bound
Male_exps2 <- Male_exps
Kegg_cf_bound2$Term <- factor(Kegg_cf_bound$Term, levels = Kegg_cf_bound.order)
kegg_rhythmic$Term <- factor(kegg_rhythmic$Term, levels = kegg_rhythmic.order)
Male_exps2$Term <- factor(Male_exps$Term, levels = Male_exps.order)
test_matrix <- Kegg_cf_bound2[,c("Term","Adjusted.P.value","Sample")]
test_matrix <- matrix(test_matrix, ncol = ncol(test_matrix), dimnames = dimnames(test_matrix)) 
#generate plot#
lengths(All_exp2)
lengths(unique(All_exp2))
test <- as.data.frame(duplicated(All_exp2$Term))
rownames(test)
test$`duplicated(All_exp2$Term)`
test_list <- lapply(test, function(x) {
 which(x == "TRUE") 
})
test_list
test_new_list <- data.frame(lapply(test_list, function(x) {
  print(All_exp2[x,])
}))
lengths(unique(test_new_list))
names(test_new_list) <- names(All_exp2)
if (x == "TRUE") {
  print() 
  } else {
  print("none")
}

lengths(All_exp)
lengths(unique(All_exp))
test <- as.data.frame(duplicated(All_exp$Term))
rownames(test)
test$`duplicated(All_exp$Term)`
test_list <- lapply(test, function(x) {
  which(x == "TRUE") 
})
test_list
test_new_list <- data.frame(lapply(test_list, function(x) {
  print(All_exp[x,])
}))
lengths(unique(test_new_list))
names(test_new_list) <- names(All_exp)
if (x == "TRUE") {
  print() 
} else {
  print("none")
}
### filters non-overlapping terms ###
test_RW_modules2 <- filter(test_RW_modules, Adjusted.P.value <= 0.05)
test_data2 <- All_exp2[duplicated(All_exp2$Term)| duplicated(All_exp2$Term, fromLast=TRUE),] 
test_data2 <- filter(test_data2, Odds.Ratio > 5)
test_data2 <- Kegg_cf_bound2[duplicated(Kegg_cf_bound2$Term)| duplicated(Kegg_cf_bound2$Term, fromLast=TRUE),] 
test_data2 <- Kegg_cm_bound2[duplicated(Kegg_cm_bound2$Term)| duplicated(Kegg_cm_bound2$Term, fromLast=TRUE),]
test_data2 <- Kegg_rm_bound2[duplicated(Kegg_rm_bound2$Term)| duplicated(Kegg_rm_bound2$Term, fromLast=TRUE),]
test_data2 <- Male_exps2[duplicated(Male_exps2$Term)| duplicated(Male_exps2$Term, fromLast=TRUE),]
test_data2 <- test_data2[duplicated(test_data2$Term)| duplicated(test_data2$Term, fromLast=TRUE),]
test_data2 <- filter(test_data2, Odds.Ratio > 4)

### clustering dot plots ###
markers <- test_data2$Term %>% unique()
markers <- Kegg_filtered$Term %>% unique()

test_data2 %>% filter(Term %in% markers) %>%
  ggplot(aes(x = Sample, y = Term, color = Adjusted.P.value, size = Odds.Ratio)) +
  geom_point() +
  scale_color_viridis_c(name = 'log2 (count +1)')

Kegg_filtered %>% filter(Term %in% markers) %>%
  ggplot(aes(x = Sample, y = Term, color = Adjusted.P.value, size = Odds.Ratio)) +
  geom_point() +
  scale_color_viridis_c(name = 'log2 (count +1)')

devtools::install_github("YuLab-SMU/ggtree")
BiocManager::install("ggdendro")
BiocManager::install("cowplot")
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork)
library(aplot)
# make data square to calculate euclidean distance
gene_cluster <- read_tsv('https://github.com/davemcg/davemcg.github.io/raw/master/content/post/scRNA_dotplot_data.tsv.gz')
gene_cluster %>% sample_n(5)
markers <- gene_cluster$Gene %>% unique()
gene_cluster %>% filter(Gene %in% markers) %>% 
  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
  geom_point() 

mat <- gene_cluster %>% 
  filter(Gene %in% markers) %>% 
  select(-cell_ct, -cell_exp_ct, -Group) %>%  # drop unused columns to faciliate widening
  pivot_wider(names_from = cluster, values_from = count) %>% 
  data.frame() # make df as tibbles -> matrix annoying
row.names(mat) <- mat$Gene  # put gene in `row`
mat <- mat[,-1] #drop gene column as now in rows
clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix


ddgram <- as.dendrogram(clust) # create dendrogram
ggtree_plot <- ggtree::ggtree(ddgram)
ggtree_plot

mat <- test_data2 %>%
  filter(Term %in% markers) %>%
  select(-Old.P.value, -Old.Adjusted.P.value, -Overlap, -P.value, -Adjusted.P.value, -Combined.Score, -Genes, -ZT, -Experiment, -lightcycle) %>% # drop unused columns to faciliate widening
  pivot_wider(names_from = Sample, values_from = Odds.Ratio) %>% 
  data.frame() # make df as tibbles -> matrix annoying
mat <- mat %>% 
  mutate_if(is.numeric, ~replace_na(., 0))
str(mat)
summary(mat)
any(is.na(mat))

row.names(mat) <- mat$Term  # put gene in `row`
mat <- mat[,-1] #drop gene column as now in rows
clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix


ddgram <- as.dendrogram(clust) # create dendrogram
ggtree_plot <- ggtree::ggtree(ddgram)
ggtree_plot

dotplot <- test_data2 %>% filter(Term %in% markers) %>%
  ggplot(aes(x = Sample, y = Term, color = Adjusted.P.value, size = Odds.Ratio)) +
  geom_point() +
  scale_color_viridis_c(name = 'log2 (count +1)') +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) #+
  #scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')

plot_grid(ggtree_plot, dotplot, nrow = 1, rel_widths = c(0.5,2), align = 'h')

dotplot2 <- test_data2 %>% 
  mutate(Gene = factor(Term, levels = clust$labels[clust$order])) %>% 
  #filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x = Sample, y = Gene, color = Adjusted.P.value, size = Odds.Ratio)) +
  geom_point() +
  cowplot::theme_cowplot() + 
  scale_color_viridis_c(name = 'log2 (count +1)') +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) #+
  #scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')

plot_grid(ggtree_plot, NULL, dotplot2, nrow = 1, rel_widths = c(0.5,-0.05, 2), align = 'h')

# make data square to calculate euclidean distance
mat <- test_data2 %>% 
  filter(Term %in% markers) %>%
  select(-Old.P.value, -Old.Adjusted.P.value, -Overlap, -P.value, -Adjusted.P.value, -Combined.Score, -Genes, -ZT, -Experiment, -lightcycle) %>% # drop unused columns to faciliate widening
  pivot_wider(names_from = Sample, values_from = Odds.Ratio) %>% 
  data.frame() # make df as tibbles -> matrix annoying
row.names(mat) <- mat$Term  # put gene in `row`
mat <- mat[,-1] #drop gene column as now in rows
mat <- mat %>% 
  mutate_if(is.numeric, ~replace_na(., 0))
v_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix
############ NOTICE THE t() above)

ddgram_col <- as.dendrogram(v_clust)
ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()

dotplot3 <- test_data2 %>% filter(Term %in% markers) %>% 
  mutate(Gene = factor(Term, levels = clust$labels[clust$order]),
         cluster = factor(Sample, levels = v_clust$labels[v_clust$order])) %>% 
  #filter(count > 0, `% Expressing` > 1) %>% 
  ggplot(aes(x = cluster, y = Gene, color = Adjusted.P.value, size = Odds.Ratio)) +
  geom_point() +
  cowplot::theme_cowplot() +
  scale_color_viridis_c(name = 'log2 (count +1)') +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  #scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)') +
  scale_y_discrete(position = "right")
#################################################
ggtree_plot_col <- ggtree_plot_col + xlim2(dotplot)

ggtree_plot <- ggtree_plot + ylim2(dotplot)
#write.csv(test_data2, "test_data2.csv")
#test_data2 <- read.csv("test_data2.csv")
library(ggplot2)
library(RColorBrewer)
nb.cols <- 24
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

labels <- ggplot(test_data2 %>% 
                   mutate(`Timepoint` = SummaryZT,
                          cluster = factor(Sample, levels = v_clust$labels[v_clust$order])), 
                 aes(x = cluster, y = 1, fill = `Timepoint`)) + 
  geom_tile() + 
  scale_fill_brewer(palette = 'Set1') + 
  #scale_fill_manual(values = mycolors) +
  theme_nothing() +
  xlim2(dotplot3)

legend <- plot_grid(get_legend(labels + theme(legend.position="bottom")))

plot_spacer() + plot_spacer() + ggtree_plot_col +
  plot_spacer() + plot_spacer() + labels + 
  plot_spacer() + plot_spacer() + plot_spacer() +
  ggtree_plot + plot_spacer() + dotplot3 + 
  plot_spacer() + plot_spacer() + legend + 
  plot_layout(ncol = 3, widths = c(0.7, -0.1, 4), heights = c(0.9, 0.1, -0.1, 4, 1))
### plot data ### 

g2 <- ggplot(kegg_nonrhythmic, aes(x = Sample, y = Term, color = Adjusted.P.value, size = Odds.Ratio)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text = element_text(face = "bold", color = "black")) +
  ylab("") +
  xlab("Non-Rhythmic Modules") +
  ggtitle("Significant KEGG Terms")
g + coord_fixed(ratio = 0.15, expand = TRUE, clip = "on")
g +  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text = element_text(face = "bold", color = "black")) + scale_y_discrete(guide = guide_axis(n.dodge = 2))
g + theme(axis.text.y = element_text(margin = margin(0,0,0,0)))
options(repr.plot.width=10, repr.plot.heigh=50)
g + theme(axis.text.y = element_text(hjust = 1))
kegg_filtered2 <- filter(Kegg_filtered, Sample != 'grey')
kegg_rhythmic <- filter(Kegg_filtered, Sample != c('grey', 'black', 'brown','darkorange', 'darkred', 'darkturquoise', 'green', 'lightcyan', 'magenta', 'purple', 'red', 'skyblue3', 'steelblue', 'turquoise', 'violet', 'yellow'))
kegg_rhythmic <- subset(Kegg_filtered, Sample == c('black', 'brown','darkorange', 'darkred', 'darkturquoise', 'green', 'lightcyan', 'magenta', 'purple', 'red', 'skyblue3', 'steelblue', 'turquoise', 'violet', 'yellow'))
kegg_rhythmic <- subset(Kegg_filtered, Sample == c("cyan", "saddlebrown", "royalblue", "skyblue", "lightyellow", "lightgreen"))
kegg_rhythmic <- Kegg_filtered %>% filter(Sample %in% c("cyan", "saddlebrown", "royalblue", "skyblue", "lightyellow", "lightgreen"))
kegg_rhythmic <- as.data.frame(kegg_rhythmic[order(kegg_rhythmic$Adjusted.P.value),])
kegg_rhythmic$Term = factor(kegg_rhythmic$Term, levels=kegg_rhythmic$Term)
kegg_nonrhythmic <- Kegg_filtered %>% filter(Sample %in% c("black", "brown","darkorange", "darkred", "darkturquoise", "green", "lightcyan", "magenta", "purple", "red", "skyblue3", "steelblue", "turquoise", "violet", "yellow"))
kegg_rhythmic.order <- unique(as.character(kegg_rhythmic$Term)[order(kegg_rhythmic$Adjusted.P.value, decreasing = TRUE)])
kegg_nonrhythmic.order <- unique(as.character(kegg_nonrhythmic$Term)[order(kegg_nonrhythmic$Adjusted.P.value, decreasing = TRUE)])
#reassign journal as factor with new levels
kegg_rhythmic$Term <- factor(kegg_rhythmic$Term, levels = kegg_rhythmic.order)
kegg_nonrhythmic$Term <- factor(kegg_nonrhythmic$Term, levels = kegg_nonrhythmic.order)
### reorder x-axis ###
g2 <- ggplot(kegg_filtered2, aes(x = factor(Sample, level = c('cyan', 'saddlebrown', 'royalblue', 'skyblue', 'lightyellow', 'lightgreen', 'black', 'brown','darkorange', 'darkred', 'darkturquoise', 'green', 'lightcyan', 'magenta', 'purple', 'red', 'skyblue3', 'steelblue', 'turquoise', 'violet', 'yellow')), y = Term, color = Adjusted.P.value, size = Odds.Ratio)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text = element_text(face = "bold", color = "black")) +
  ylab("") +
  xlab("Modules") +
  ggtitle("Significant KEGG Terms")
## reorder x-axis using factor function ## 
ggplot(test_data2, aes(x = factor(Sample, level = c('CM_salmon', 'CM_purple', 'RW_blue', 'RW_cyan', 'RW_magenta', 'CM_royalblue', 'CM_magenta', 'RW_red', 'CM_brown','RW_brown', 'RW_pink', 'RW_purple', 'CM_cyan', 'RW_lightcyan', 'RW_royalblue', 'CM_midnightblue', 'RW_green', 'CM_green', 'RW_yellow', 'CM_lightcyan', 'CM_red', 'RW_darkred', 'CM_yellow', 'RW_black', 'RW_midnightblue', 'RW_tan', 'RW_lightgreen', 'RW_darkgrey')), y = Term, color = Adjusted.P.value, size = Odds.Ratio)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text = element_text(face = "bold", color = "black")) +
  ylab("") +
  xlab("Modules") +
  ggtitle("Significant KEGG Terms")

library(plotly)
ggplotly(p)
###sep###
names(test_kg) <- modules_interest
names(test_kg)
setwd("/Users/aron/Desktop/LaSalle_Lab/Analysis/clamsrw/rnaseq/separate_normalization/GOTermsUpdatedCorrected")

data = read.csv(glue::glue("yellow_moduleRM12.csv")) 
RM_yellow = read.csv(glue::glue("yellow_moduleRM12.csv")) 
RM_black = read.csv(glue::glue("black_moduleRM12.csv"))
RM_pink = read.csv(glue::glue("pink_moduleRM12.csv"))
RM_red = read.csv(glue::glue("red_moduleRM12.csv"))
Red_RW <- RM_red %>%
    dplyr::select(x) %>%
    purrr::flatten() %>%
    enrichR::enrichr(c("GO_Biological_Process_2023",
                       "GO_Cellular_Component_2023",
                       "GO_Molecular_Function_2023",
                       "KEGG_2019_Mouse",
                       "Panther_2016",
                       "Reactome_2022",
                       "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO")) %>% 
    purrr::set_names(names(.) %>% stringr::str_trunc(40, ellipsis = "")) %T>%
    purrr::map(~ dplyr::filter(., Adjusted.P.value < 0.05))
#%>% 
    #purrr::map(~ dplyr::filter(., stringr::str_detect(Genes, ";")))
    #openxlsx::write.xlsx(file = glue::glue("Module_{module}_RM12_enrichr.xlsx")) 

RW_Yellow_BP <- Yellow_RW$GO_Biological_Process_2023
RW_Yellow_BP$Sample <- "RW_yellow"
RW_Black_BP <- Black_RW$GO_Biological_Process_2023 
RW_Black_BP$Sample <- "RW_black"
test_RW_modules <- rbind(RW_Yellow_BP, RW_Black_BP)
head(test_RW_modules)
test_RW_modules2 <- filter(test_RW_modules, Adjusted.P.value <= 0.05)

RW_Yellow_KG <- Yellow_RW$KEGG_2019_Mouse
RW_Yellow_KG$Sample <- "RW_yellow"
RW_Pink_KG <- Pink_RW$KEGG_2019_Mouse
RW_Pink_KG$Sample <- "RW_pink"
RW_Red_KG <- Red_RW$KEGG_2019_Mouse
RW_Red_KG$Sample <- "RW_red"
test_RW_modules <- rbind(RW_Yellow_KG, RW_Pink_KG, RW_Red_KG)
head(test_RW_modules)
test_RW_modules2 <- filter(test_RW_modules, Adjusted.P.value <= 0.05)
#gobp1 <- test9$GO_Biological_Process_2023[test9$GO_Biological_Process_2023$Adjusted.P.value < 0.05,]
ggplot(test_RW_modules2, aes(x = Sample, y = Term, color = Adjusted.P.value, size = Odds.Ratio)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() +
  ylab("") +
  xlab("Modules") +
  ggtitle("GO Terms")


####
    DMRichR::slimGO(tool = "enrichR",
                    annoDb = "org.Mm.eg.db",
                    plots = FALSE) %T>%
    openxlsx::write.xlsx(file = glue::glue("Module_{module}_WT121_rrvgo_enrichr.xlsx")) %>%
    DMRichR::GOplot() %>%
    ggplot2::ggsave(glue::glue("Module_{module}_WT121_enrichr_plot.pdf"),
                    plot = .,
                    device = NULL,
                    height = 8.5,
                    width = 10) 
  


test2 <- read.csv("module_blue.csv")
