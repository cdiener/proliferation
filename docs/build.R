#!/usr/bin/env Rscript

setwd("..")
rmarkdown::render("protocol.rmd", output_format="html_document", output_file="docs/index.html")
gc()
rmarkdown::render("protocol.rmd", output_format="pdf_document", output_file="docs/protocol.pdf")
setwd("docs")
