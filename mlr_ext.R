#  mlr_ext.R
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

# custom regressor for the H2O deep learner
makeRLearner.regr.h2odeep = function() {
    makeRLearnerRegr(
    cl = "regr.h2odeep",
    package = "h2o",
    par.set = makeParamSet(
        makeDiscreteLearnerParam(id="activation", default="RectifierWithDropout",
        values=c("Rectifier", "Tanh", "TanhWithDropout", "RectifierWithDropout", 
        "Maxout", "MaxoutWithDropout")),
        makeIntegerVectorLearnerParam(id="hidden", default=c(200, 200)),
        makeNumericLearnerParam(id="epochs", default=10, lower=0),
        makeNumericLearnerParam(id="rate", lower=0),
        makeNumericLearnerParam(id="l1", lower=0),
        makeNumericLearnerParam(id="l2", lower=0)
    ),
    properties = c("numerics", "factors"),
    name = "H2O Deep Learning",
    short.name = "h2odeep",
    note = ""
    )
}

trainLearner.regr.h2odeep = function(.learner, .task, .subset, .weights = NULL, ...) {
    h2o.init(startH2O=F)
    dat <- getTaskData(.task, .subset)
    features <- names(dat)
    target <- .task$task.desc$target
    features <- features[features != target]
    tf <- as.h2o(dat)

    h2o::h2o.deeplearning(x=features, y=target, training_frame=tf, 
        quiet_mode=T, shuffle_training_data=T, ...)
}

predictLearner.regr.h2odeep = function(.learner, .model, .newdata, ...) {
  nd <- as.h2o(.newdata)
  as.data.frame(h2o.predict(.model$learner.model, newdata = nd))$predict
}

# custom measure for median relative error
mre.fun = function(task, model, pred, feats, extra.args) {
  ae = abs(getPredictionResponse(pred) - getPredictionTruth(pred))
  median(ae/abs(getPredictionTruth(pred)))
}

mre = makeMeasure(
  id = "mre", name = "Mean relative error",
  properties = c("regr", "req.pred", "req.truth"),
  minimize = TRUE, best = 0, worst = Inf,
  fun = mre.fun
)

# Helper function to train and predict in a single step
train_perf <- function(learner, task) {
    tr <- train(learner, task)
    p <- predict(tr, task)
    list(model=tr, pred=p, perf=performance(p, measures=list(mse, mae, 
        mre, rsq)))
}

combine_measures <- function(results, extract="measures.test", 
    method=c("elastic net", "random forest", "gradient boosting", 
    "deep learning")) {
    n <- length(method)
    per_met <- lapply(1:n, function(i) cbind(results[[i]][[extract]], i))
    summ <- data.table::rbindlist(per_met)
    summ$i <- method[summ$i]
    names(summ)[ncol(summ)] <- "method"
    summ
}
