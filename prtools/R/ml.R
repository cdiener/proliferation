# Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
# MIT license. See LICENSE for more information.

#' Various measures for goodness-of-fit.
#'
#' Calculates a set of five different goodness-of-fit metrics.
#'
#' @param truth A numeric vector describing the true values.
#' @param pred A numeric vector describing the predicted values.
#' @return A numeric vector of length 5 containing the following metrics:
#' \describe{
#'  \item{mse}{The mean squared error.}
#'  \item{rmse}{The root mean squared error.}
#'  \item{mae}{The mean absolute error.}
#'  \item{mre}{The mean relative error.}
#'  \item{rsq}{The r-squared value.}
#' }
#' @examples
#'  NULL
#'
#' @export
measures <- function(truth, pred) {
    res <- truth - pred
    mse <- mean(res ^ 2)
    rmse <- sqrt(mse)
    mae <- mean(abs(res))
    mre <- mean(abs(res) / abs(truth))
    rsq <- 1 - sum(res ^ 2) / sum( (truth - mean(truth) ) ^ 2)

    c(mse = mse, rmse = rmse, mae = mae, mre = mre, rsq = rsq)
}

#' Calculates interaction terms for a gene expression matrix.
#'
#' Interaction terms for two expression values g[i] and g[j] are given
#' by g[i]*g[j].
#'
#' @param M A gene expression matrix with genes in the columns and samples
#'  in the rows.
#' @return A new gene expression matrix with n + 0.5*n*(n-1) columns, where n
#'  is the number of columns in the matrix M.
#' @examples
#'  NULL
#'
#' @export
inter <- function(M) {
    nam <- colnames(M)
    out <- NULL
    for (i in 1:ncol(M)) {
        newnames <- paste0(nam[i], "x", colnames(M[, i:ncol(M), drop = FALSE]))
        x <- M[, i] * M[, i:ncol(M), drop = FALSE]
        colnames(x) <- newnames
        out <- cbind(out, x)
    }

    return(out)
}
