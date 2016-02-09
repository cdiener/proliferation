#  tcga.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

devtools::load_all("~/code/tcgar")
library(data.table)
library(ggplot2)

folders <- dir("TCGA", pattern="[A-Z]{2,4}$", include.dirs=T)
folders <- file.path("TCGA", folders)

tcga <- read_bulk(folders)
gc()

all_rnaseq <- rowMeans(log(tcga$RNASeqV2$counts+1))
all_huex <- rowMeans(tcga$HuEx$assay)

nci60 <- fread("regprob.csv", header=T)
rates <- nci60$rates
nci60[, rates := NULL]
all_nci60 <- colMeans(nci60)
imp <- fread("rfimportance.csv")
imp <- imp[importance > sqrt(.Machine$double.eps)]
imp <- imp[importance > quantile(importance, .99)]
imp <- sub("hsa_", "", imp$entrez)

cat("Comparing HuEx\n")
data(huex_bm)
setkey(huex_bm, entrez)
shared <- intersect(names(all_nci60), huex_bm$entrez)
huex <- data.frame(NCI60=all_nci60[huex_bm[shared, entrez]], TCGA=all_huex[huex_bm[shared, symbol]])
#huex_imp <- data.frame(NCI60=all_nci60[imp$entrez], TCGA=all_huex[huex_features[imp$entrez, symbol]])

huex_plot <- ggplot(huex, aes(x=TCGA, y=NCI60)) + geom_point(alpha=0.25) + 
    geom_abline() #+ geom_point(data=huex_imp, col="red")

cat("Comparing HuEx")
shared <- intersect(names(all_nci60), names(all_rnaseq))
rna <- data.frame(NCI60=all_nci60[shared], TCGA=all_rnaseq[shared])
#rna_imp <- data.frame(NCI60=all_nci60[imp$entrez], TCGA=all_rnaseq[imp$entrez])

rna_plot <- ggplot(rna, aes(x=TCGA, y=NCI60)) + geom_point(alpha=0.25) + 
    geom_abline()
