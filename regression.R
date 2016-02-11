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

measures = function(truth, pred) {
    res <- truth - pred
    mse <- mean(res^2)
    rmse <- sqrt(mse)
    mae <- mean(abs(res))
    mre <- mean(abs(res)/abs(truth))
    rsq <- 1 - sum(res^2)/sum((truth - mean(truth))^2)
    
    c(mse=mse, rmse=rmse, mae=mae, mre=mre, rsq=rsq)
}

inter = function(M) {
    nam <- colnames(M)
    out <- NULL
    for(i in 1:ncol(M)) {
        newnames <- paste0(nam[i], "x", colnames(M[, i:ncol(M), drop=F]))
        x <- M[, i]*M[, i:ncol(M), drop=F]
        colnames(x) <- newnames
        out <- cbind(out, x)
    }
    
    return(out)
}

start_t <- proc.time()
registerDoParallel(cl=6)

cat("Reading data\n")
rdata <- fread("regprob.csv", header=T)
load("rna_calib.rda")
good <- rna$ensgene[(rna$TCGA - rna$NCI60)^2 < 4]
#names(rdata)[-ncol(rdata)] <- paste0("hsa_", names(rdata)[-ncol(rdata)])
rates <- rdata$rates
rdata[, "rates" := NULL]
rdata <- as.matrix(rdata)[, good]

cat("Running 1st order regressor\n")

pred <- data.frame()
m <- data.frame()

mod1 <- cv.glmnet(rdata, rates, nfolds=length(rates), keep=T, parallel=T)
pred_train <- predict(mod1, rdata, s="lambda.min")[,1]
pred_test <- mod1$fit.preval[, which.min(mod1$cvm)]
pred <- rbind(pred, data.frame(truth=rates, pred=pred_train, set="train", order="1st"))
pred <- rbind(pred, data.frame(truth=rates, pred=pred_test, set="test", order="1st"))
m <- rbind(m, data.frame(t(measures(rates, pred_train)), set="train", order="1st"))
m <- rbind(m, data.frame(t(measures(rates, pred_test)), set="test", order="1st"))
nonzero <- abs(coef(mod1, s="lambda.min")[-1]) > 0

cat("Running 2nd order regressor\n")

data2 <- inter(rdata[,nonzero])
mod2 <- cv.glmnet(data2, rates, nfolds=10, keep=T, parallel=T)
pred_train <- predict(mod2, data2, s="lambda.min")[,1]
pred_test <- mod2$fit.preval[, which.min(mod2$cvm)]
pred <- rbind(pred, data.frame(truth=rates, pred=pred_train, set="train", order="2nd"))
pred <- rbind(pred, data.frame(truth=rates, pred=pred_test, set="test", order="2nd"))
m <- rbind(m, data.frame(t(measures(rates, pred_train)), set="train", order="2nd"))
m <- rbind(m, data.frame(t(measures(rates, pred_test)), set="test", order="2nd"))

cat("Running 1st + 2nd order regressor\n")

data12 <- cbind(rdata[, nonzero], data2)
mod3 <- cv.glmnet(data12, rates, nfolds=length(rates), keep=T, parallel=T)
pred_train <- predict(mod3, data12, s="lambda.min")[,1]
pred_test <- mod3$fit.preval[, which.min(mod3$cvm)]
pred <- rbind(pred, data.frame(truth=rates, pred=pred_train, set="train", order="1st and 2nd"))
pred <- rbind(pred, data.frame(truth=rates, pred=pred_test, set="test", order="1st and 2nd"))
m <- rbind(m, data.frame(t(measures(rates, pred_train)), set="train", order="1st and 2nd"))
m <- rbind(m, data.frame(t(measures(rates, pred_test)), set="test", order="1st and 2nd"))

cat("Reducing model by cutoff\n")
cf <- as.numeric(coef(mod2, s="lambda.min"))[-1]
names(cf) <- rownames(coef(mod2))[-1]
nonzero <- abs(cf) > 1e-3
data_red <- data2[, nonzero]
mod <- cv.glmnet(data_red, rates, nfolds=length(rates), keep=T, parallel=T)
pred_train <- predict(mod, data_red, s="lambda.min")[,1]
pred_test <- mod$fit.preval[, which.min(mod$cvm)]
pred <- rbind(pred, data.frame(truth=rates, pred=pred_train, set="train", order="2nd cutoff"))
pred <- rbind(pred, data.frame(truth=rates, pred=pred_test, set="test", order="2nd cutoff"))
m <- rbind(m, data.frame(t(measures(rates, pred_train)), set="train", order="2nd cutoff"))
m <- rbind(m, data.frame(t(measures(rates, pred_test)), set="test", order="2nd cutoff"))

genes <- do.call(rbind, strsplit(colnames(data_red), "x"))
colnames(genes) <- c("gene1", "gene2")
readr::write_csv(data.frame(genes, coef=cf[nonzero]), "best_interactions.csv")
save(mod, file="glmnet_model.rda")

#Assemble predictions
pred_plot <- ggplot(pred, aes(x=truth, y=pred, col=order)) + geom_abline() + 
    geom_point() + facet_grid(set ~ order) + theme_bw()

write("----------\nUsed time:\n----------", file = "")
print(proc.time() - start_t) 
