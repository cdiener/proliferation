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
    start <- proc.time()
    tcga <- read_bulk(folders)
    write("------Read data------", file="")
    print(proc.time() - start)
    gc()
    save(tcga, file="tcga.rda")
} else load("tcga.rda")

all_rnaseq <- log(rowMeans(tcga$RNASeqV2$counts)+1, 2)
all_rnaseq <- data.table(entrez=tcga$RNASeqV2$features$entrez,
    tcga_rnaseq=all_rnaseq)
all_huex <- data.table(symbol=tcga$HuEx$features$symbol,
    tcga_huex=rowMeans(tcga$HuEx$assay))

nci60 <- fread("regprob.csv", header=T)
rates <- nci60$rates
nci60[, rates := NULL]
all_nci60 <- data.table(ensgene=colnames(nci60), nci60_huex=colMeans(nci60))

ints <- fread("best_interactions.csv", header=T)
genes <- unique(c(ints$gene1, ints$gene2))

join <- merge(as.data.table(rnaseq_bm), all_nci60, by="ensgene")
join <- merge(join, all_huex, by="symbol")
join <- merge(join, all_rnaseq, by="entrez")
dupes <- join$ensgene[duplicated(join$ensgene)]
join <- join[!(join$ensgene %in% dupes), ]
imp <- join[ensgene %in% genes]

cat("Comparing HuEx - NCI60\n")
norm_int <- mean(join$nci60_huex) - mean(join$tcga_huex)

huex_plot <- ggplot(join, aes(x=tcga_huex, y=nci60_huex)) + geom_point(alpha=0.1, stroke=0) +
    geom_abline(intercept=norm_int, size=3/4, col="dodgerblue") +
    geom_abline(intercept=sqrt(2)+norm_int, size=3/4, col="dodgerblue", linetype="dashed") +
    geom_abline(intercept=-sqrt(2)+norm_int, size=3/4, col="dodgerblue", linetype="dashed") +
    geom_point(data=imp, col="red", stroke=0) + xlab("TCGA HuEx") +
    ylab("NCI60 HuEx") + theme_bw()
ggsave("images/huex.png", huex_plot, width=100, height=90, units="mm", dpi=300)

cat("Comparing RNASeq - HuEx\n")
norm <- coef(glm(tcga_huex ~ tcga_rnaseq + 0, data=join))

intern_plot <- ggplot(join, aes(x=tcga_rnaseq, y=tcga_huex)) + geom_point(alpha=0.1, stroke=0) +
    geom_abline(slope=norm, size=3/4, col="dodgerblue") +
    geom_point(data=imp, col="red", stroke=0) + xlab("TCGA RNA-seq") +
    ylab("TCGA HuEx") + theme_bw()
ggsave("images/tcga.png", intern_plot, width=100, height=90, units="mm", dpi=300)

cat("Comparing RNASeq - NCI60\n")

rna_plot <- ggplot(join, aes(x=tcga_rnaseq, y=nci60_huex)) + geom_point(alpha=0.1, stroke=0) +
    geom_abline(slope=norm, intercept=norm_int, size=3/4, col="dodgerblue") +
    geom_abline(slope=norm, intercept=sqrt(2)+norm_int, size=3/4, col="dodgerblue", linetype="dashed") +
    geom_abline(slope=norm, intercept=-sqrt(2)+norm_int, size=3/4, col="dodgerblue", linetype="dashed") +
    geom_point(data=imp, col="red", stroke=0) + xlab("TCGA RNA-seq") +
    ylab("NCI60 HuEx") + theme_bw()
ggsave("images/rnaseq.png", rna_plot, width=100, height=90, units="mm", dpi=300)

cat("Extracting joined samples...\n")
save(join, file="join.rda")
save(norm, norm_int, file="norm_factor.rda")
