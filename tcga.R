#  tcga.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

devtools::load_all("~/code/tcgar")
library(data.table)
library(ggplot2)

cat("Reading TCGA data..\n")

folders <- dir("TCGA", pattern="[A-Z]{2,4}$", include.dirs=T)
folders <- file.path("TCGA", folders)

if (!file.exists("tcga.rda")) {
    tcga <- read_bulk(folders)
    gc()
    save(tcga, file="tcga.rda")
} else load("tcga.rda")

all_rnaseq <- rowMeans(log(tcga$RNASeqV2$counts+1))
all_huex <- rowMeans(tcga$HuEx$assay)

nci60 <- fread("regprob.csv", header=T)
rates <- nci60$rates
nci60[, rates := NULL]
all_nci60 <- colMeans(nci60)
ints <- fread("best_interactions.csv", header=T)
genes <- unique(c(ints$gene1, ints$gene2))

cat("Comparing HuEx\n")
hdf <- data.table(huex_bm)
setkey(hdt, ensgene)
map <- hdt[shared, list(ensgene, symbol)]
shared <- intersect(names(all_nci60), hdt$ensgene)
huex <- data.frame(NCI60=all_nci60[map$ensgene], TCGA=all_huex[map$symbol])
save(huex, file="huex_calib.rda")
map <- hdt[genes, list(ensgene, symbol)]
huex_imp <- data.frame(NCI60=all_nci60[map$ensgene], TCGA=all_huex[map$symbol], map)

huex_plot <- ggplot(huex, aes(x=TCGA, y=NCI60)) + geom_point(alpha=0.25) + 
    geom_abline() + geom_point(data=huex_imp, col="red")

cat("Comparing RNASeq\n")
rdt <- data.table(rnaseq_bm)
setkey(rdt, ensgene)
shared <- intersect(names(all_nci60), rdt$ensgene)
map <- rdt[shared, list(ensgene, entrez)]
rna <- data.frame(NCI60=all_nci60[map$ensgene], TCGA=all_rnaseq[map$entrez], map)
save(rna, file="rna_calib.rda")
map <- rdt[genes, list(ensgene, entrez)]
rna_imp <- data.frame(NCI60=all_nci60[map$ensgene], TCGA=all_rnaseq[map$entrez])

rna_plot <- ggplot(rna, aes(x=TCGA, y=NCI60)) + geom_point(alpha=0.25) + 
    geom_abline()+ geom_point(data=rna_imp, col="red")
