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

# For bootstrapping or CV
stats <- function(x, i, alpha, lambda, y) {
    fit <- glmnet(x[i,], y[i], alpha=alpha)
    ip <- if (length(y[-i])>0) -i else i
    pred <- predict(fit, x[ip,], s=lambda)[,1]
    mse <- mean((y[ip] - pred)^2)
    rel_error <- mean(abs(y[ip] - pred)/abs(0.5*(y[ip]+pred)))
    co <- cor(y[ip], pred)
    return(c(mse=mse, rel_error=rel_error, cor=co))
}

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
    i <- 0
    n <- length(unique(groups))
    ex <- exprs(eset)
    env <- environment()
    
    summ <- function(ids) {
        if (progress) {
            cat("                                                           \r")
            cat(sprintf("Joining group %d/%d...", env$i <- env$i + 1, n))
        }
        fun(ex[ids,,drop=F])
    }
    
    res <- tapply(features, groups, summ, simplify=FALSE)
    res <- do.call(rbind, res)
    cat("\n")
    
    return(ExpressionSet(res))
}

### RUN ###

cat("Loading ExpressionSet and assigning annotations...\n")
if (file.exists("eset.Rd")) load("eset.Rd") else {
    celfiles <- list.celfiles(recursive=T, listGzipped=T)
    raw_data <- read.celfiles(celfiles)
    rm(raw_data)
    eset <- rma(raw_data, target="probeset")
    genemap <- fread("huex_to_enstxp.csv", colClasses=rep("character", 3))
    setkey(genemap, huex)
    
    eset <- eset[rownames(eset) %in% genemap$huex, ]
    genemap <- genemap[huex %in% rownames(eset)]

    cat("Reducing ExpressionSet...\n")
    # Summarize by transcript
    eset <- eset_reduce(eset, genemap$huex, genemap$enstxp, colMeans)
    save(eset, file="eset.Rd")
}

samples <- fread("samples.csv")

cat("Preparing data for growth rate inference...\n")
# Correlate growth rates to transcripts
gcs <- fread("growth_rates.csv")
setkey(gcs, cell_line)
cell_lines <- intersect(gcs$cell_line, samples$cell_line)
eset_summ <- sapply(cell_lines, function(cl) rowMeans(exprs(eset)[,samples$cell_line==cl]))
colnames(eset_summ) <- cell_lines
rates <- gcs[cell_lines, log(2)/doubling_time*24]
names(rates) <- cell_lines
if (!all(colnames(eset_summ) == names(rates))) cat("Wrong cell line ordering!")

# cors <- apply(eset_summ, 1, cor, y=rates)

# Do elastic net regression
cat("Getting elastic net GLM...\n")
library(glmnet)

# Grid search for good alpha strategy
al <- seq(0.1,1,by=0.05)
mods <- lapply(al, function(a) cv.glmnet(t(eset_summ), log(rates), alpha=a, 
    parallel=TRUE, nfolds=10))
res <- lapply(mods, function(m) m$cvm)
x <- NULL
for (i in 1:length(al)) {
    x <- rbind(x, data.frame(mse=res[[i]], alpha=al[i]))
}

alpha_plot <- ggplot(x, aes(x=alpha, y=mse)) + geom_jitter() + stat_smooth()
best <- which.min(sapply(res, min))
cat(paste0("Found best model at alpha=", al[best], "\n"))
mod <- mods[[best]]

pdf("alpha_steps.pdf", width=8, height=6)
plot(mod)
dev.off()

pred <- predict(mod, t(eset_summ), s=1e-2)

df <- data.frame(predicted=exp(pred[,1]), measured=rates, panel=gcs[cell_lines, panel_name])
fit_plot <- ggplot(df, aes(x=predicted, y=measured)) + 
    geom_point(aes(col=panel)) + geom_abline()




 
