## Copyright 2016 Christian Diener <mail[at]cdiener.com>
##
## MIT license. See LICENSE for more information.

tissue_lfc <- function(v, map, extra) {
    p <- map[rownames(v)]

    out <- NULL

    out <- lapply(unique(p), function(pname) {
        ti <- p == pname
        in_t <- colMeans(v[ti,]) + 1e-6
        out_t <- colMeans(v[!ti,]) + 1e-6
        lfcs <- log(in_t, 2) - log(out_t, 2)

        res <- data.table(rid=names(in_t), panel=pname, lfc=lfcs)
        cbind(res, extra)
    })

    return(rbindlist(out))
}

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

shorten <- function(x, n=32) {
    sapply(x, function(xi) {
        xi <- as.character(xi)
        if (nchar(xi) > n) {
            xi <- paste0(substr(xi,1,n-1), "...")
        }
        xi
    })
}
