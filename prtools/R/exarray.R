## Copyright 2016 Christian Diener <mail[at]cdiener.com>
##
## MIT license. See LICENSE for more information.

#' @useDynLib prtools
#' @importFrom Rcpp sourceCpp
NULL

#' Reduces an ExpressionSet by an n-to-n map of features to groups.
#'
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
#'
#' @export
#' @importFrom Biobase ExpressionSet exprs
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
