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

nfr <- mean(all_rnaseq)/mean(all_nci60)

cat("Comparing HuEx - NCI60\n")
hdt <- data.table(huex_bm)
setkey(hdt, ensgene)
shared <- intersect(names(all_nci60), hdt$ensgene)
map <- hdt[shared, list(ensgene, symbol)]
huex <- data.frame(NCI60=all_nci60[map$ensgene], TCGA=all_huex[map$symbol], map)
map <- hdt[genes, list(ensgene, symbol)]
huex_imp <- data.frame(NCI60=all_nci60[map$ensgene], TCGA=all_huex[map$symbol], map)

huex_plot <- ggplot(huex, aes(x=TCGA, y=NCI60)) + geom_point(alpha=0.25) + 
    geom_smooth(method="glm", formula=y ~ x + 0) + 
    geom_point(data=huex_imp, col="red") + theme_bw()
ggsave("huex.png", huex_plot, width=6, height=6, dpi=300)

cat("Comparing RNASeq - NCI60\n")
rdt <- data.table(rnaseq_bm)
setkey(rdt, ensgene)
shared <- intersect(names(all_nci60), rdt$ensgene)
map <- rdt[shared, list(ensgene, entrez)]
rna <- data.frame(NCI60=all_nci60[map$ensgene], TCGA=all_rnaseq[map$entrez], map)
map <- rdt[genes, list(ensgene, entrez)]
rna_imp <- data.frame(NCI60=all_nci60[map$ensgene], TCGA=all_rnaseq[map$entrez], map)

rna_plot <- ggplot(rna, aes(x=TCGA, y=NCI60)) + geom_point(alpha=0.25) + 
    geom_smooth(method="glm", formula=y ~ x + 0) + 
    geom_point(data=rna_imp, col="red") + theme_bw()
ggsave("rnaseq.png", rna_plot, width=6, height=6, dpi=300)

cat("Comparing RNASeq - HuEx\n")
setkey(hdt, entrez)
shared <- intersect(hdt$entrez, rdt$entrez)
map <- hdt[shared, list(entrez, symbol)]
intern <- data.frame(huex=all_huex[map$symbol], rnaseq=all_rnaseq[map$entrez], map)
setkey(hdt, ensgene)
map <- hdt[genes, list(entrez, symbol)]
intern_imp <- data.frame(huex=all_huex[map$symbol], rnaseq=all_rnaseq[map$entrez], map)

intern_plot <- ggplot(intern, aes(x=rnaseq, y=huex)) + geom_point(alpha=0.25) + 
    geom_smooth(method="glm", formula=y ~ x + 0) + 
    geom_point(data=intern_imp, col="red") + theme_bw()
ggsave("TCGA_control.png", intern_plot, width=6, height=6, dpi=300)

cat("Extracting joined samples...\n")
join <- merge(rna, huex, by=c("ensgene", "NCI60"), suffixes=c("_RNA", "_MA"))
dupes <- join$ensgene[duplicated(join$ensgene)]
join <- join[!(join$ensgene %in% dupes), ]
save(join, file="join.rda")

cat("Getting normalization factors...\n")
norm_nci <- coef(glm(TCGA ~ NCI60 + 0, data=rna))
norm_control <- coef(glm(rnaseq ~ huex + 0, data=intern))

save(norm_nci, file="norm_factor.rda")
