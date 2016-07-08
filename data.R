## Copyright 2016 Christian Diener <mail[at]cdiener.com>
##
## MIT license. See LICENSE for more information.

library(prtools)

# Download the NCI-60 HuEx data
GEOquery::getGEOSuppFiles("GSE29682")
untar("GSE29682/GSE29682_RAW.tar", exdir="GSE29682")

# Get the TCGA data
create.dir("TCGA")
setwd("TCGA")
panels <- tcgar::get_panels()
none <- sapply(panels, tcgar::get_data, tech=c("clinical", "HuEx", "RNASeqV2"))
setwd("..")

# Get annotations for genes
library(biomaRt)
library(data.table)

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
attrs <- c("ensembl_gene_id", "description", "external_gene_name", "ucsc",
    "entrezgene", "affy_huex_1_0_st_v2")
genemap <- getBM(mart=ensembl, attributes=attrs)
genemap <- data.table(genemap)
names(genemap) <- c("ensgene", "description", "symbol", "ucsc", "entrez", "huex")
saveRDS(genemap, "genemap.rds")
