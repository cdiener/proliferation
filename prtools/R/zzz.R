# Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
# MIT license. See LICENSE for more information.

#' @import RcppProgress
#' @importFrom utils packageVersion
NULL

pkgs <- c("data.table", "oligo", "foreach", "doMC", "survival", "tcgar",
    "glmnet", "biomaRt", "ggplot2", "pheatmap")

silent_lib <- function(...) suppressPackageStartupMessages(library(...))

.onAttach <- function(...) {
    is_loaded <- paste0("package:", pkgs) %in% search()
    needed <- sort(pkgs[!is_loaded])

    if (length(needed) == 0) return()

    vs <- sapply(needed, function(x) as.character(packageVersion(x)))
    packageStartupMessage(paste("Also loading:", needed, vs, collapse = "\n"))
    lapply(needed, silent_lib, character.only = TRUE, warn.conflicts = FALSE)
}
