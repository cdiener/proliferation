#  tcga.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

library(prtools)

# Read TCGA data
tcga <- list(
    rnaseq = read_rnaseq("GDC/manifest_tcga_rnaseq.tsv", "GDC/rnaseq"),
    huex = read_huex("GDC/manifest_tcga_huex.tsv", "GDC/huex"),
    clinical = read_clinical("GDC/manifest_tcga_clinical.tsv", "GDC/clinical")
)
saveRDS(tcga, "tcga.rds")

all_rnaseq <- log(rowMeans(tcga$rnaseq$counts) + 1, 2)
all_rnaseq <- data.table(ensgene = tcga$rnaseq$features$ensgene,
                     tcga_rnaseq = all_rnaseq)
all_huex <- data.table(symbol = tcga$huex$features$symbol,
                    tcga_huex = rowMeans(tcga$huex$assay))

nci60 <- fread("regprob.csv", header = T)
rates <- nci60$rates
nci60[, rates := NULL]
all_nci60 <- data.table(ensgene = colnames(nci60), nci60_huex = colMeans(nci60))

join <- merge(as.data.table(genemap), all_nci60, by = "ensgene")
join <- merge(join, all_huex, by = "symbol")
join <- merge(join, all_rnaseq, by = "ensgene")
dupes <- join$ensgene[duplicated(join$ensgene)]
join <- join[!(join$ensgene %in% dupes), ]

ints <- fread("best_interactions.csv", header = T)
genes <- unique(c(ints$gene1, ints$gene2))
imp <- join[ensgene %in% genes]

cat("Comparing HuEx - NCI60\n")
norm_int <- mean(join$nci60_huex) - mean(join$tcga_huex)

huex_plot <- ggplot(join, aes(x = tcga_huex, y = nci60_huex)) +
    geom_point(alpha = 0.1, stroke = 0) +
    geom_abline(intercept = norm_int, size = 0.75, col = "dodgerblue") +
    geom_abline(intercept = sqrt(2) + norm_int, size = 0.75, col = "dodgerblue",
    linetype = "dashed") + geom_abline(intercept = -sqrt(2) + norm_int,
    size = 0.75, col = "dodgerblue", linetype = "dashed") +
    geom_point(data = imp, col = "red", stroke = 0) + xlab("TCGA HuEx") +
    ylab("NCI60 HuEx") + theme_bw()
ggsave("images/huex.png", huex_plot, width = 100, height = 90, units = "mm",
    dpi = 300)

cat("Comparing RNASeq - HuEx\n")
norm <- coef(glm(tcga_huex ~ tcga_rnaseq + 0, data = join))

intern_plot <- ggplot(join, aes(x = tcga_rnaseq, y = tcga_huex)) +
    geom_point(alpha = 0.1, stroke = 0) +
    geom_abline(slope = norm, size = 0.75, col = "dodgerblue") +
    geom_point(data = imp, col = "red", stroke = 0) + xlab("TCGA RNA-seq") +
    ylab("TCGA HuEx") + theme_bw()
ggsave("images/tcga.png", intern_plot, width = 100, height = 90, units = "mm",
    dpi = 300)

cat("Comparing RNASeq - NCI60\n")

rna_plot <- ggplot(join, aes(x = tcga_rnaseq, y = nci60_huex)) +
    geom_point(alpha = 0.1, stroke = 0) +
    geom_abline(slope = norm, intercept = norm_int, size = 0.75,
        col = "dodgerblue") +
    geom_abline(slope = norm, intercept = sqrt(2) + norm_int, size = 0.75,
        col = "dodgerblue", linetype = "dashed") +
    geom_abline(slope = norm, intercept = -sqrt(2) + norm_int, size = 0.75,
        col = "dodgerblue", linetype = "dashed") +
    geom_point(data = imp, col = "red", stroke = 0) + xlab("TCGA RNA-seq") +
    ylab("NCI60 HuEx") + theme_bw()
ggsave("images/rnaseq.png", rna_plot, width = 100, height = 90, units = "mm",
    dpi = 300)

cat("Extracting joined samples...\n")
saveRDS(join, "join.rds")
saveRDS(c(norm, norm_int), file = "norm_factor.rds")
