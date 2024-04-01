###2.20.24
cf <- read.csv("cf.counts.csv")
cm <- read.csv("cm.counts.csv")
rf <- read.csv("rf.counts.csv")
rm <- read.csv("rm.counts.csv")
wtmf <- read.csv("wtfm.csv")
wtf <- read.csv("wtf_non-normalized.csv")
wtm <- read.csv("wtm_non-normalized.csv")
wtf.1 <- read.csv("wt.fcounts.csv")
wtm.1 <- read.csv("wt.mcounts.csv")
####

rownames(cf) <- cf$Gene_ID
cf <- cf[,-c(1)]

rownames(cm) <- cm$Gene_ID
cm <- cm[,-c(1)]

rownames(rf) <- rf$Gene_ID
rf <- rf[,-c(1)]

rownames(rm) <- rm$Gene_ID
rm <- rm[,-c(1)]

rownames(wtf) <- wtf$Gene_ID
wtf <- wtf[, -c(1)]

rownames(wtm) <- wtm$Gene_ID
wtm <- wtm[, -c(1)]

### Set up DGEList
d0 <- DGEList(cf)
cutoff <- 2
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(cf) # number of genes before cleanup
dim(d) # number of genes left
logcpm <- cpm(d, prior.count=2, log=TRUE)
write.csv(logcpm,"cf_normalized.csv")
cf_normalized <- logcpm

d0 <- DGEList(cm)
cutoff <- 2
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(cm) # number of genes before cleanup
dim(d) # number of genes left
logcpm <- cpm(d, prior.count=2, log=TRUE)
write.csv(logcpm,"cm_normalized.csv")
cm_normalized <- logcpm

d0 <- DGEList(rf)
cutoff <- 2
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(rf) # number of genes before cleanup
dim(d) # number of genes left
logcpm <- cpm(d, prior.count=2, log=TRUE)
write.csv(logcpm,"rf_normalized.csv")
rf_normalized <- logcpm

d0 <- DGEList(rm)
cutoff <- 2
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(rm) # number of genes before cleanup
dim(d) # number of genes left
logcpm <- cpm(d, prior.count=2, log=TRUE)
write.csv(logcpm,"rm_normalized.csv")
rm_normalized <- logcpm

d0 <- DGEList(wtf)
cutoff <- 2
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(wtf) # number of genes before cleanup
dim(d) # number of genes left
logcpm <- cpm(d, prior.count=2, log=TRUE)
write.csv(logcpm,"wtf_normalized.csv")
wtf_normalized <- logcpm
wtf_normalized <- wtf_normalized[, order(colnames(wtf_normalized))]

d0 <- DGEList(wtm)
cutoff <- 2
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(wtm) # number of genes before cleanup
dim(d) # number of genes left
logcpm <- cpm(d, prior.count=2, log=TRUE)
write.csv(logcpm,"wtm_normalized.csv")
wtm_normalized <- logcpm
wtm_normalized <- wtm_normalized[, order(colnames(wtm_normalized))]

###
wtfwtm <- merge(wtf_normalized, wtm_normalized, by = "row.names")
colnames(wtfwtm)[1] <- "Gene_ID"
cfcm <- merge(cf_normalized, cm_normalized, by = "row.names")
rfrm <- merge(rf_normalized, rm_normalized, by = "row.names")
cfcmrfrmwtfm <- merge(wtfwtm, cfcmrfrm_normalized, by = "Gene_ID")
cfcmrfrm_normalized <- merge(cfcm, rfrm, by = "Row.names")
cfcmrfrmwtfm_normalized <- merge(wtmf, cfcmrfrm_normalized, by = "Row.names")
All_merge <- merge(wtfwtm, cfcmrfrm_normalized, by = "Gene_ID")
dim(cfcmrfrmwtfm_normalized)
write.csv(cfcmrfrmwtfm_normalized, "cfcmrfrmwtfm_normalized.csv")

### separate out data
rownames(All_merge) <- All_merge$Gene_ID
All_merge <- All_merge[,-c(1)]
head(All_merge)
wtf_norm <- All_merge[,c(1:46)]
wtm_norm <- All_merge[,c(47:92)]
cf_norm <- All_merge[,c(93:124)]
cm_norm <- All_merge[,c(125:155)]
rf_norm <- All_merge[,c(156:187)]
rm_norm <- All_merge[,c(188:217)]
colnames(rm_norm)

### Merge Traits
CLAMS_alltraits <- read.csv("CLAMS_alltraits.csv")
RER_traits <- read.csv("CLAMS 12-11 LD-DD combined traits.csv")
merge_RERtraits <- merge(CLAMS_alltraits, RER_traits, by = "ID")
