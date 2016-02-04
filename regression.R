#  regression.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

library(data.table)
library(mlr)
library(randomForestSRC)
library(parallelMap)
library(h2o)

source("mlr_ext.R")

start_t <- proc.time()

cat("Reading data\n")
rdata <- fread("regprob.csv", header=T)
names(rdata)[-ncol(rdata)] <- paste0("hsa_", names(rdata)[-ncol(rdata)])

cat("Running full regressors\n")
task <- makeRegrTask(id = "nci60", data=as.data.frame(rdata), target="rates")
task <- normalizeFeatures(task, method="standardize")
lrn <- list(makeLearner("regr.glmnet", par.vals=list(alpha=0.5)), 
    makeLearner("regr.randomForestSRC", par.vals=list(ntree=50)), 
    makeLearner("regr.xgboost", par.vals=list(max_depth=2, eta=1e-1, 
    nrounds=100, verbose=0)),
    makeLearner("regr.h2odeep", par.vals=list(hidden=c(100, 100, 100), epochs=10,
    l1=0.01, l2=0.01)))
resa <- makeResampleDesc(method = "CV", iters=8)

parallelStartSocket(7, show.info=F)
parallelLibrary("mlr", "h2o")
parallelExport("makeRLearner.regr.h2odeep", "trainLearner.regr.h2odeep",
    "predictLearner.regr.h2odeep", "mre")
h2o.init(nthreads=7)

r <- benchmark(learner = lrn, task = task, resampling = resa, 
    measures=list(mse, mae, mre), show.info=F)
train_full <- parallelMap(train_perf, lrn, more.args=list(task=task))


cat("Getting feature importance\n")
rfimp <- rfsrc(rates ~ ., data=rdata, ntree=1000, importance="permute.ensemble")$importance
cutoff <- quantile(rfimp[rfimp > sqrt(.Machine$double.eps)], .99)
imp <- names(rdata)[rfimp > cutoff]
cat(sprintf("- kept %d variables\n", length(imp)))

cat("Running reduced regressors\n")
redtask <- subsetTask(task, features=imp)

redlrn <- list(makeLearner("regr.glmnet", par.vals=list(alpha=0.5)), 
    makeLearner("regr.randomForestSRC", par.vals=list(ntree=1000, mtry=length(imp))), 
    makeLearner("regr.xgboost", par.vals=list(max_depth=2, eta=1e-2, 
    nrounds=1e2, verbose=0)),
    makeLearner("regr.h2odeep", par.vals=list(hidden=c(100, 100, 100), 
    l1=0.01, epochs=10)))
redsa <- makeResampleDesc(method = "LOO")

r_red <- benchmark(learner = redlrn, task = redtask, resampling = redsa, 
    measures=list(mse, mae, mre), show.info=F)
train_red <- parallelMap(train_perf, redlrn, more.args=list(task=redtask))

#Assemble predictions
preds_test_full <- combine_measures(r$results$nci60, extract=c("pred", "data"))
test_full_plot <- ggplot(preds_test_full, aes(x=truth, y=response, col=method)) + 
    geom_abline() + geom_point() + facet_wrap(~method, nrow=1) + theme_bw()

preds <- combine_measures(train_red, extract=c("pred", "data"))
train_plot <- ggplot(preds, aes(x=truth, y=response, col=method)) + geom_abline() + 
    geom_point() + facet_wrap(~method, nrow=1) + theme_bw()

measures <- combine_measures(r_red$results$nci60)
measures <- melt(measures, id.vars=c("iter", "method"))
cv_plot <- ggplot(measures, aes(x=method, y=value, col=method)) + geom_boxplot() +
    geom_jitter() + facet_wrap(~variable, scales="free") + theme_bw()
    
preds_test <- combine_measures(r_red$results$nci60, extract=c("pred", "data"))
preds_test <- rbindlist(preds_test)
preds_test$i <- method[preds_test$i]
names(preds_test)[4] <- "method"
test_plot <- ggplot(preds_test, aes(x=truth, y=response, col=method)) + geom_abline() + 
    geom_point() + facet_wrap(~method, nrow=1) + theme_bw()

parallelStop()


write("----------\nUsed time:\n----------", file = "")
print(proc.time() - start_t) 
