#  regression.R
#
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#
#  MIT license. See LICENSE for more information.

library(data.table)
library(doParallel)
library(foreach)
library(glmnet)
library(ggplot2)
library(tcgar)
library(prtools)

start_t <- proc.time()
registerDoParallel(cl=6)

cat("Reading data\n")
rdata <- fread("regprob.csv", header=T)
load("norm_factor.rda")
load("join.rda")
good <- join[(tcga_rnaseq*norm + norm_int - nci60_huex)^2 < 1 &
    (tcga_huex + norm_int - join$nci60_huex)^2 < 1, ensgene]
#names(rdata)[-ncol(rdata)] <- paste0("hsa_", names(rdata)[-ncol(rdata)])
rates <- rdata$rates
rdata[, "rates" := NULL]
rdata <- as.matrix(rdata)[, good]

cat("Running 1st order regressor\n")

pred <- data.frame()
m <- data.frame()

mod1 <- cv.glmnet(rdata, rates, nfolds=length(rates), keep=T, parallel=T,
    grouped=FALSE, standardize=FALSE)
pred_train <- predict(mod1, rdata, s="lambda.min")[,1]
pred_test <- mod1$fit.preval[, which.min(mod1$cvm)]
pred <- rbind(pred, data.frame(truth=rates, pred=pred_train, set="train", order="1st"))
pred <- rbind(pred, data.frame(truth=rates, pred=pred_test, set="validation", order="1st"))
m <- rbind(m, data.frame(t(measures(rates, pred_train)), set="train", order="1st"))
m <- rbind(m, data.frame(t(measures(rates, pred_test)), set="validation", order="1st"))
nonzero <- abs(coef(mod1, s="lambda.min")[-1]) > 0

cat("Running 2nd order regressor\n")

data2 <- inter(rdata[,nonzero])
mod2 <- cv.glmnet(data2, rates, nfolds=10, keep=T, parallel=T,
    grouped=FALSE, standardize=FALSE)
pred_train <- predict(mod2, data2, s="lambda.min")[,1]
pred_test <- mod2$fit.preval[, which.min(mod2$cvm)]
pred <- rbind(pred, data.frame(truth=rates, pred=pred_train, set="train", order="2nd"))
pred <- rbind(pred, data.frame(truth=rates, pred=pred_test, set="validation", order="2nd"))
m <- rbind(m, data.frame(t(measures(rates, pred_train)), set="train", order="2nd"))
m <- rbind(m, data.frame(t(measures(rates, pred_test)), set="validation", order="2nd"))

cat("Running 1st + 2nd order regressor\n")

data12 <- cbind(rdata[, nonzero], data2)
mod3 <- cv.glmnet(data12, rates, nfolds=length(rates), keep=T, parallel=T,
    grouped=FALSE, standardize=FALSE)
pred_train <- predict(mod3, data12, s="lambda.min")[,1]
pred_test <- mod3$fit.preval[, which.min(mod3$cvm)]
pred <- rbind(pred, data.frame(truth=rates, pred=pred_train, set="train", order="1st and 2nd"))
pred <- rbind(pred, data.frame(truth=rates, pred=pred_test, set="validation", order="1st and 2nd"))
m <- rbind(m, data.frame(t(measures(rates, pred_train)), set="train", order="1st and 2nd"))
m <- rbind(m, data.frame(t(measures(rates, pred_test)), set="validation", order="1st and 2nd"))

cat("Reducing model by cutoff\n")
cf <- as.numeric(coef(mod2, s="lambda.min"))[-1]
names(cf) <- rownames(coef(mod2))[-1]
nonzero <- abs(cf) > 1e-5
data_red <- data2[, nonzero]
mod <- cv.glmnet(data_red, rates, nfolds=length(rates), keep=T,
    grouped=FALSE, standardize=FALSE)
pred_train <- predict(mod, data_red, s="lambda.min")[,1]
pred_test <- mod$fit.preval[, which.min(mod$cvm)]
pred <- rbind(pred, data.frame(truth=rates, pred=pred_train, set="train", order="2nd + cutoff"))
pred <- rbind(pred, data.frame(truth=rates, pred=pred_test, set="validation", order="2nd + cutoff"))
m <- rbind(m, data.frame(t(measures(rates, pred_train)), set="train", order="2nd + cutoff"))
m <- rbind(m, data.frame(t(measures(rates, pred_test)), set="validation", order="2nd + cutoff"))

genes <- do.call(rbind, strsplit(colnames(data_red), "x"))
colnames(genes) <- c("gene1", "gene2")
readr::write_csv(data.frame(genes, coef=cf[nonzero]), "best_interactions.csv")
save(mod, file="glmnet_model.rda")

#Assemble predictions
pred_plot <- ggplot(pred, aes(x=truth, y=pred, col=order)) + geom_abline() +
    geom_point() + facet_grid(set ~ order) + theme_bw() +
    xlab("measured proliferation rate [1/h]") +
    ylab("predicted proliferation rate [1/h]") +
    theme(legend.position="none")
ggsave("images/model.png", pred_plot, width=185, height=90, units="mm", dpi=300)

cat("Predicting...\n")
load("tcga.rda")
hdt <- unique(data.table(huex_bm), by="ensgene")
setkey(hdt, ensgene)

symbs <- cbind(hdt[genes[,1], symbol], hdt[genes[,2], symbol])
huex_ex <- tcga$HuEx$assay + norm_int
huex_red <- t(huex_ex[symbs[,1], ]*huex_ex[symbs[,2], ])
colnames(huex_red) <- paste0(symbs[,1], "x", symbs[,2])
rates_huex <- predict(mod, huex_red, s="lambda.min")[,1]
controls <- is.na(names(rates_huex))
rates_huex <- rates_huex[!controls]

rdt <- unique(data.table(rnaseq_bm), by="ensgene")
setkey(rdt, ensgene)

entrez <- cbind(rdt[genes[,1], entrez], rdt[genes[,2], entrez])
rna_ex <- log(tcga$RNASeqV2$counts+1, 2)
rna_ex <- rna_ex * norm + norm_int
rna_red <- t(rna_ex[entrez[,1], ]*rna_ex[entrez[,2], ])
colnames(rna_red) <- paste0(entrez[,1], "x", entrez[,2])
rates_rna <- predict(mod, rna_red, s="lambda.min")[,1]
tum <- tcga$RNASeqV2$samples$tumor

pred <- data.table(barcode=c(names(rates_rna), names(rates_huex)),
    rates=c(rates_rna, rates_huex),
    panel=c(tcga$RNASeqV2$samples$panel, tcga$HuEx$samples$panel[!controls]),
    tumor=c(tcga$RNASeqV2$samples$tumor, tcga$HuEx$samples$tumor[!controls]))
#pred <- unique(pred, by="barcode")

comb <- merge(pred, tcga$clinical$samples, by=c("barcode", "tumor", "panel"))
comb <- merge(comb, tcga$clinical$patient, by=c("patient_uuid", "patient_barcode"))

save(pred, comb, file="combined.rda")
readr::write_csv(pred, "pred_rates.csv")

write("----------\nUsed time:\n----------", file = "")
print(proc.time() - start_t)
