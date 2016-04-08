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

all_rnaseq <- log(rowMeans(tcga$RNASeqV2$counts)+1, 2)
all_rnaseq <- all_rnaseq
all_huex <- rowMeans(tcga$HuEx$assay)

nci60 <- fread("regprob.csv", header=T)
rates <- nci60$rates
nci60[, rates := NULL]
all_nci60 <- colMeans(nci60)
ints <- fread("best_interactions.csv", header=T)
genes <- unique(c(ints$gene1, ints$gene2))

cat("Comparing HuEx - NCI60\n")
hdt <- data.table(huex_bm)
setkey(hdt, ensgene)
shared <- intersect(names(all_nci60), hdt$ensgene)
map <- hdt[shared, list(ensgene, symbol)]
huex <- data.frame(NCI60=all_nci60[map$ensgene], TCGA=all_huex[map$symbol], map)
map <- hdt[genes, list(ensgene, symbol)]
huex_imp <- data.frame(NCI60=all_nci60[map$ensgene], TCGA=all_huex[map$symbol], map)
norm_int <- mean(huex$NCI60) - mean(huex$TCGA)

huex_plot <- ggplot(huex, aes(x=TCGA, y=NCI60)) + geom_point(alpha=0.1, stroke=0) +
    geom_abline(intercept=norm_int, size=3/4, col="dodgerblue") +
    geom_abline(intercept=sqrt(2)+norm_int, size=3/4, col="dodgerblue", linetype="dashed") +
    geom_abline(intercept=-sqrt(2)+norm_int, size=3/4, col="dodgerblue", linetype="dashed") +
    geom_point(data=huex_imp, col="red", stroke=0) + xlab("TCGA HuEx") +
    ylab("NCI60 HuEx") + theme_bw()
ggsave("images/huex.png", huex_plot, width=100, height=90, units="mm", dpi=300)

cat("Comparing RNASeq - HuEx\n")
setkey(hdt, entrez)
shared <- intersect(hdt$entrez, rdt$entrez)
map <- hdt[shared, list(entrez, symbol)]
intern <- data.frame(huex=all_huex[map$symbol], rnaseq=all_rnaseq[map$entrez], map)
setkey(hdt, ensgene)
map <- hdt[genes, list(entrez, symbol)]
norm <- coef(glm(huex ~ rnaseq + 0, data=intern))
intern_imp <- data.frame(huex=all_huex[map$symbol], rnaseq=all_rnaseq[map$entrez], map)

intern_plot <- ggplot(intern, aes(x=rnaseq, y=huex)) + geom_point(alpha=0.1, stroke=0) +
    geom_abline(slope=norm, size=3/4, col="dodgerblue") +
    geom_point(data=intern_imp, col="red", stroke=0) + xlab("TCGA RNA-seq") +
    ylab("TCGA HuEx") + theme_bw()
ggsave("images/tcga.png", intern_plot, width=100, height=90, units="mm", dpi=300)

cat("Comparing RNASeq - NCI60\n")
rdt <- data.table(rnaseq_bm)
setkey(rdt, ensgene)
shared <- intersect(names(all_nci60), rdt$ensgene)
map <- rdt[shared, list(ensgene, entrez)]
rna <- data.frame(NCI60=all_nci60[map$ensgene], TCGA=all_rnaseq[map$entrez], map)
map <- rdt[genes, list(ensgene, entrez)]
rna_imp <- data.frame(NCI60=all_nci60[map$ensgene], TCGA=all_rnaseq[map$entrez], map)

rna_plot <- ggplot(rna, aes(x=TCGA, y=NCI60)) + geom_point(alpha=0.1, stroke=0) +
    geom_abline(slope=norm, intercept=norm_int, size=3/4, col="dodgerblue") +
    geom_abline(slope=norm, intercept=sqrt(2)+norm_int, size=3/4, col="dodgerblue", linetype="dashed") +
    geom_abline(slope=norm, intercept=-sqrt(2)+norm_int, size=3/4, col="dodgerblue", linetype="dashed") +
    geom_point(data=rna_imp, col="red", stroke=0) + xlab("TCGA RNA-seq") +
    ylab("NCI60 HuEx") + theme_bw()
ggsave("images/rnaseq.png", rna_plot, width=100, height=90, units="mm", dpi=300)

cat("Extracting joined samples...\n")
join <- merge(rna, huex, by=c("ensgene", "NCI60"), suffixes=c("_RNA", "_MA"))
dupes <- join$ensgene[duplicated(join$ensgene)]
join <- join[!(join$ensgene %in% dupes), ]
save(join, file="join.rda")

save(norm, norm_int, file="norm_factor.rda")
