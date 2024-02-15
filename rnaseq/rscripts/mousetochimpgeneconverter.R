setwd("/Users/aron/Downloads/HumanChimp-Data (3)")
getwd()
file = bzfile("/Users/aron/Downloads/HumanChimp-Data (3)/Dataset 1 (network construction).csv.bz2")
dat1 = read.csv(file, header = T)

datExpr=data.frame(t(dat1[dat1$Brain_variant_H>0,2:39]))
indexHuman=c(19:36)
indexChimp=c(1:18)
nSets=2
multiExpr=list()
multiExpr[[1]]=list(data=datExpr[indexHuman,])
multiExpr[[2]]=list(data=datExpr[indexChimp,]) 
colnames(multiExpr[[1]]$data)=dat1$Probe_set[dat1$Brain_variant_H>0] 
colnames(multiExpr[[2]]$data)=dat1$Probe_set[dat1$Brain_variant_H>0] 
probes=colnames(multiExpr[[1]]$data)

BiocManager::install("KEGG.db")

BiocManager::install("hgu95av2.db")

checkSets(exp)$nGenes

nGenes = checkSets(exp)$nGenes
setNames = c("wt", "crw")
mp = list()
cr = list()
doModulePreservation = TRUE
doClusterRepro = TRUE
ensemblCodes = colnames(exp$femdata$data)
head(ensemblCodes)
anno <- read.delim("ensembl_mm_100.tsv",as.is=T)
mapped_probes <- mappedkeys(anno)
data2db = match(ensemblCodes, anno$Gene.stable.ID)
fin = is.finite(data2db)
head(data2db)
test <- as.data.frame(data2db)

library(hgu95av2.db)
x<-hgu95av2ENTREZID
mapped_probes<-mappedkeys(x)
xx <- as.list(x[mapped_probes])
dbProbes=names(xx); 
dbEntrez=sapply(xx,I)
data2db=match(probes,dbProbes); 
fin=is.finite(data2db); 
probeEntrez=dbEntrez[data2db[fin]] 
multiExpr[[1]]$data=multiExpr[[1]]$data[,fin]; 
multiExpr[[2]]$data=multiExpr[[2]]$data[,fin];
x<-hgu95av2SYMBOL
mapped_probes<-mappedkeys(x)
xx<-as.list(x[mapped_probes]) 
dbProbes=names(xx); 
dbGeneNames=sapply(xx,I) 
data2db=match(probes,dbProbes); 
fin=is.finite(data2db); 
probeGeneNames=dbGeneNames[data2db[fin]] 
entrez2name = cbind(probeEntrez,probeGeneNames);

colMedians=function(x){apply(x,2,median,na.rm=TRUE)} 
geneExpr=list(); 
for(set in 1:nSets) 
{ 
  x=collapseRows(t(multiExpr[[set]]$data),probeEntrez,probes[fin], 
                 method="function",methodFunction=colMedians, 
                 connectivityPower=9); 
  geneExpr[[set]]=list(data=t(x$datETcollapsed)); 
}
