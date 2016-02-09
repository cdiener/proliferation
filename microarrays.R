#  microarrays.R
#  
#  Copyright 2015 Christian Diener <ch.diener[a]gmail.com>
#  
#  MIT license. See LICENSE for more information.

library(foreach)
library(doParallel)
registerDoParallel(cl=6)
library(oligo)
library(stringr)
library(data.table)
library(ggplot2)

#' Reduces an ExpressionSet by an n-to-n map of features to groups. All entries
#' in \code{features} must exist in \code{eset}. \code{features} and 
#' \code{groups} must have the same length.
#'
#' @param eset An ExpressionSet object.
#' @param features A character vector of features to be grouped.
#' @param groups A factor or character vector mapping the entries in 
#'  \code{features} to groups.
#' @param fun The reduction function. Must operate on a n x m matrix and return
#'  a vector of length m.
#' @param progress Should progress information be shown.
#' @return A new ExpressionSet with features given by \code{unique(groups)}.
eset_reduce <- function(eset, features, groups, fun, progress=TRUE) {
    groups <- factor(groups)
    ex <- exprs(eset)
    
    idx <- 1:nrow(ex)
    names(idx) <- rownames(ex)
    idx <- idx[features]
    facs <- as.integer(groups)
    res <- matrix_reduce(ex, idx, facs, "mean")
    rownames(res) <- levels(groups)
    
    return(ExpressionSet(res))
}

### RUN ###

cat("Loading ExpressionSet and assigning annotations...\n")
if (file.exists("eset_raw.Rd")) load("eset_raw.Rd") else {
    celfiles <- list.celfiles(recursive=T, listGzipped=T)
    raw_data <- read.celfiles(celfiles)
    eset <- rma(raw_data, target="probeset")
    rm(raw_data)
    save(eset, "eset_raw.Rd")
}

genemap <- fread("genemap.csv", colClasses=rep("character", 7), 
    na.strings=c("NA",""))
genemap <- unique(genemap, by=c("huex", "entrez"))
setkey(genemap, huex)

eset <- eset[rownames(eset) %in% genemap$huex, ]
genemap <- genemap[huex %in% rownames(eset) & !is.na(entrez)]
Rcpp::sourceCpp("matrix_reduce.cpp")
cat("Reducing ExpressionSet...\n")
# Summarize by Ensembl Gene
eset <- eset_reduce(eset, genemap$huex, genemap$entrez)
save(eset, file="eset.Rd")

samples <- fread("samples.csv")

cat("Preparing data for growth rate inference...\n")
gcs <- fread("growth_rates.csv")
setkey(gcs, cell_line)
cell_lines <- intersect(gcs$cell_line, samples$cell_line)
eset_summ <- sapply(cell_lines, function(cl) rowMeans(exprs(eset)[,samples$cell_line==cl]))
colnames(eset_summ) <- cell_lines
rates <- gcs[cell_lines, log(2)/doubling_time*24]
names(rates) <- cell_lines
if (!all(colnames(eset_summ) == names(rates))) cat("Wrong cell line ordering!")

# cors <- apply(eset_summ, 1, cor, y=rates)

# use only features that have shared Ensembl Gene IDs in NCI60 and TCGA
devtools::load_all("~/code/tcgar")
data(huex_features)
data(rnaseqv2_features)
#shared <- intersect(rownames(eset_summ), rnaseqv2_features$entrez)

export <- data.table(t(eset_summ))
export[, "rates" := rates]
readr::write_csv(export, "regprob.csv")
#source("regression.R")




 
