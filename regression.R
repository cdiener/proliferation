#  regression.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

library(data.table)
library(mlr)
library(randomForestSRC)
library(parallelMap)

train_perf <- function(learner, task) {
    p <- predict(train(learner, task), task)
    performance(p, measures=list(mse, medae, rsq))
}

start_t <- proc.time()

cat("Reading data\n")
rdata <- fread("regprob.csv", header=T)
names(rdata)[-ncol(rdata)] <- paste0("hsa_", names(rdata)[-ncol(rdata)])

cat("Running full regressors\n")
task <- makeRegrTask(id = "nci60", data=as.data.frame(rdata), target="rates")
task <- normalizeFeatures(task, method="standardize")
lrn <- list(makeLearner("regr.glmnet", par.vals=list(lambda=0.2)), 
    makeLearner("regr.randomForestSRC", par.vals=list(ntree=50)), 
    makeLearner("regr.xgboost", par.vals=list(max_depth=2, eta=1e-1, 
    nrounds=100, verbose=0)))
resa <- makeResampleDesc(method = "CV")

parallelStartSocket(7, show.info=F)
parallelLibrary("mlr")

r <- benchmark(learner = lrn, task = task, resampling = resa, 
    measures=list(mse, medae), show.info=F)
train_full <- parallelMap(train_perf, lrn, more.args=list(task=task))

cat("Getting feature importance\n")
rfimp <- rfsrc(rates ~ ., data=rdata, ntree=5000)$importance
cutoff <- quantile(rfimp[rfimp > sqrt(.Machine$double.eps)], .99)
imp <- names(rdata)[rfimp > cutoff]
cat(sprintf("- kept %d variables\n", length(imp)))

cat("Running reduced regressors\n")
redtask <- subsetTask(task, features=imp)

redlrn <- list(makeLearner("regr.glmnet", par.vals=list(lambda=0.2)), 
    makeLearner("regr.randomForestSRC", par.vals=list(ntree=1000, mtry=length(imp))), 
    makeLearner("regr.xgboost", par.vals=list(max_depth=2, eta=1e-2, 
    nrounds=1e3, verbose=0)))
redsa <- makeResampleDesc(method = "RepCV", reps=10)

r_red <- benchmark(learner = lrn, task = redtask, resampling = redsa, 
    measures=list(mse, medae), show.info=F)
train_red <- parallelMap(train_perf, redlrn, more.args=list(task=redtask))

parallelStop()

write("----------\nUsed time:\n----------", file = "")
print(proc.time() - start_t) 
