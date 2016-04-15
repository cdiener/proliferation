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

#' Reduces an ExpressionSet by an n-to-n map of features to groups.

#' All entries in \code{features} must exist in \code{eset}. \code{features} and
#' \code{groups} must have the same length. Note that the 'mean' and 'median'
#' are geometric variants meaning they operate on the log-expression values
#' rather than on the raw expression values which usually gives better results.
#'
#' @param eset An ExpressionSet object.
#' @param features A character vector of features to be grouped.
#' @param groups A factor or character vector mapping the entries in
#'  \code{features} to groups.
#' @param mean The reduction method. Must be one of 'mean', 'median' or 'max'.
#' @param progress Should progress information be shown.
#' @return A new ExpressionSet with features given by \code{unique(groups)}.
eset_reduce <- function(eset, features, groups, method="mean", progress=TRUE) {
    groups <- factor(groups)
    ex <- exprs(eset)

    if (!(method %in% c("mean", "median", "max"))) stop("Not a valid method!")

    idx <- 1:nrow(ex)
    names(idx) <- rownames(ex)
    idx <- idx[features]
    facs <- as.integer(groups)

    res <- matrix_reduce(ex, idx, facs, method)
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
    save(eset, file="eset_raw.rda")
}

genemap <- fread("genemap.csv", colClasses=rep("character", 7),
    na.strings=c("NA",""))
genemap <- unique(genemap, by=c("huex", "ensgene"))
setkey(genemap, huex)

eset <- eset[rownames(eset) %in% genemap$huex, ]
genemap <- genemap[huex %in% rownames(eset) & !is.na(ensgene)]
Rcpp::sourceCpp("matrix_reduce.cpp")
cat("Reducing ExpressionSet...\n")
# Summarize by Ensembl Gene
eset <- eset_reduce(eset, genemap$huex, genemap$ensgene)
save(eset, file="eset.rda")

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
