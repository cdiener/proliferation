## Copyright 2016 Christian Diener <mail[at]cdiener.com>
##
## MIT license. See LICENSE for more information.

#' Calculates flux log-fold changes between in-tissue and out-tissue samples.
#'
#' `tissue_lfc` takes values for fluxes across different tissues and calculates
#' the log-fold changes between the flux in a given target tissue against all
#' all other tissues, for all tissues.
#'
#' @param v A matrix of fluxes where rows denote samples across different
#' tissues and columns denote fluxes.
#' @param map A map translating the row names of v to their specific tissue.
#' @param extra Optional data.frame with as many rows as rows in v that contains
#' additional information about the fluxes that should be appended to the results.
#' @return A data.table containing the log-fold changes for each flux and tissue.
#' @examples
#'  NULL
#'
#' @export
#' @importFrom data.table data.table rbindlist
tissue_lfc <- function(v, map, extra=NULL) {
    p <- map[rownames(v)]

    out <- NULL

    out <- lapply(unique(p), function(pname) {
        ti <- p == pname
        in_t <- colMeans(v[ti,]) + 1e-6
        out_t <- colMeans(v[!ti, ]) + 1e-6
        lfcs <- log(in_t, 2) - log(out_t, 2)

        res <- data.table(rid=names(in_t), panel=pname, lfc=lfcs)
        if (!is.null(extra)) res <- cbind(res, extra)
        res
    })

    return(rbindlist(out))
}

#' Calculates GSEA-like enrichment scores.
#'
#' @param p The name of the pathway for which the enrichment score is to be
#' calculated.
#' @param w A ranked vector of values/weights.
#' @param pws A vector with as many entries as in w, assigning a pathway to each
#' entry in w.
#' @param both Optional logical parameter. Should the enrichment score be
#' separated by negative/positive values.
#' @return The enrichment score.
#' @examples
#'  NULL
#'
#' @export
ES <- function(p, w, pws, both=FALSE) {
    n <- length(pws)
    nr <- sum(abs(w[pws == p]))
    nh <- sum(pws == p)

    scores <- vector(length=n)
    scores[pws == p] <- abs(w[pws == p])/nr
    scores[pws != p] <- -1/(n - nh)
    r <- range(cumsum(scores))
    i <- which.max(abs(r))
    if (both) i <- 1:2

    r[i]
}

#' Calculates the normalized enrichment scores.
#'
#' Calculates es/mean(perm), where perm is a set of n random permutations
#' of pathway labels.
#'
#' @param p The name of the pathway for which the enrichment score is to be
#' calculated.
#' @param w A ranked vector of values/weights.
#' @param pws A vector with as many entries as in w, assigning a pathway to each
#' entry in w.
#' @param n The number of random permutations to generate.
#' @return A vector with two elements: the normalized enrichment score and an
#' empirical, uncorrected p-value for it.
#' @examples
#'  NULL
#'
#' @export
NES <- function(p, w, pws, n=200) {
    es <- ES(p, w, pws)
    norm <- replicate(n, ES(p, w, sample(pws), both=TRUE))

    if (es < 0) {
        no <- norm[norm < 0]
        nes <- -es/mean(no)
        pval <- sum(no < es)/length(no)
    } else {
        no <- norm[norm >= 0]
        nes <- es/mean(no)
        pval <- sum(no > es)/length(no)
    }

    return(c(nes, pval))
}

#' Shorten a string to given length.
#'
#' @param x A string or vector of strings to be shortened.
#' @param n The maximum length of the resulting string(s).
#' @return The shortened string(s).
#' @examples
#' print(shorten("bla"))
#' print(shorten("blablablablablablablablablablablablablablablabla"))
#'
#' @export
shorten <- function(x, n=32) {
    short <- sapply(x, function(xi) {
        xi <- as.character(xi)
        if (nchar(xi) > n) {
            xi <- paste0(substr(xi,1,n-1), "...")
        }
        xi
    })

    unname(short)
}
