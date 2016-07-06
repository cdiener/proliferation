#  microarrays.R
#
#  Copyright 2015 Christian Diener <ch.diener[a]gmail.com>
#
#  MIT license. See LICENSE for more information.

library(prtools)
library(oligo)
library(data.table)
library(foreach)

doParallel::registerDoParallel(cl=6)

cat("Loading ExpressionSet and assigning annotations...\n")
if (file.exists("eset_raw.rds")) eset_raw <- readRDS("eset_raw.rda") else {
    celfiles <- list.celfiles(recursive=T, listGzipped=T)
    raw_data <- read.celfiles(celfiles)
    eset <- rma(raw_data, target="probeset")
    rm(raw_data)
    saveRDS(eset, file="eset_raw.rds")
}

if (file.exists("genemap.rda")) genemap <- readRDS("genemap.rds") else {
    ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
    attrs <- c("ensembl_gene_id", "description", "external_gene_name", "ucsc",
        "entrezgene", "affy_huex_1_0_st_v2")
    genemap <- getBM(mart=ensembl, attributes=attrs)
    genemap <- data.table(genemap)
    names(genemap) <- c("ensgene", "description", "symbol", "ucsc", "entrez", "huex")
    saveRDS(genemap, "genemap.rds")
}

genemap <- unique(genemap, by=c("ensgene", "huex"))
setkey(genemap, huex)
eset <- eset[rownames(eset) %in% genemap$huex, ]
genemap <- genemap[huex %in% rownames(eset) & !is.na(ensgene)]

cat("Reducing ExpressionSet...\n")
# Summarize by Ensembl Gene
eset <- eset_reduce(eset, genemap$huex, genemap$ensgene)

samples <- fread("samples.csv")

cat("Preparing data for growth rate inference...\n")
gcs <- fread("growth_rates.csv")
setkey(gcs, cell_line)
cell_lines <- intersect(gcs$cell_line, samples$cell_line)
eset_summ <- sapply(cell_lines, function(cl) rowMeans(exprs(eset)[,samples$cell_line==cl]))
colnames(eset_summ) <- cell_lines
rates <- gcs[cell_lines, log(2)/doubling_time]
names(rates) <- cell_lines
if (!all(colnames(eset_summ) == names(rates))) cat("Wrong cell line ordering!")

# cors <- apply(eset_summ, 1, cor, y=rates)


export <- data.table(t(eset_summ))
export[, "rates" := rates]
readr::write_csv(export, "regprob.csv")
#source("regression.R")
