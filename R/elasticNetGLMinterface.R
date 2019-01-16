setGeneric("elasticNetGLMtrainInterface", function(measurements, ...)
{standardGeneric("elasticNetGLMtrainInterface")})

setMethod("elasticNetGLMtrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  elasticNetGLMtrainInterface(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

# Clinical data or one of the other inputs, transformed.
setMethod("elasticNetGLMtrainInterface", "DataFrame", function(measurements, classes, lambda = NULL, ..., verbose = 3)
{
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting elastic net regularised GLM classifier to data.")
  
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- data.frame(splitDataset[["measurements"]], check.names = FALSE)
  measurementsMatrix <- model.matrix(splitDataset[["classes"]] ~ -1 + ., data = measurements, check.names = FALSE)
  
  fitted <- glmnet::glmnet(measurementsMatrix, splitDataset[["classes"]], family = "multinomial", type.multinomial = "grouped", ...)

  if(is.null(lambda)) # fitted has numerous models for automatically chosen lambda values.
  { # Pick one lambda based on resubstitution performance.
    bestLambda <- fitted[["lambda"]][which.min(sapply(fitted[["lambda"]], function(lambda) # Largest Lambda with minimum balanced error rate.
    {
      classPredictions <- factor(as.character(predict(fitted, measurementsMatrix, s = lambda, type = "class")), levels = fitted[["classnames"]])
      calcExternalPerformance(classes, classPredictions, "balanced error")
    }))[1]]
    attr(fitted, "tune") <- list(lambda = bestLambda)
  }

  attr(fitted, "features") <- colnames(measurements)[attr(measurementsMatrix, "assign")]
  fitted
})

# One or more omics datasets, possibly with clinical data.
setMethod("elasticNetGLMtrainInterface", "MultiAssayExperiment",
function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  
  if(ncol(measurements) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    elasticNetGLMtrainInterface(measurements, classes, ...)
})

# Matrix of numeric measurements.
setGeneric("elasticNetGLMpredictInterface", function(model, test, ...)
{standardGeneric("elasticNetGLMpredictInterface")})

setMethod("elasticNetGLMpredictInterface", c("multnet", "matrix"),
          function(model, test, ...)
{
  elasticNetGLMpredictInterface(model, DataFrame(t(test), check.names = FALSE), ...)
})

# Clinical data only.
setMethod("elasticNetGLMpredictInterface", c("multnet", "DataFrame"), function(model, test, classes = NULL, lambda, ..., returnType = c("class", "score", "both"), verbose = 3)
{ # ... just consumes emitted tuning variables from .doTrain which are unused.
  if(!is.null(classes))
  {
    splitDataset <- .splitDataAndClasses(test, classes)  # Remove any classes, if present.
    test <- splitDataset[["measurements"]]
  }
  
  returnType <- match.arg(returnType)
  
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained elastic net regularised GLM classifier.")

  if(missing(lambda)) # Tuning parameters are not passed to prediction functions.
    lambda <- attr(model, "tune")[["lambda"]] # Sneak it in as an attribute on the model.

  classPredictions <- factor(as.character(predict(model, as.matrix(test), s = lambda, type = "class")), levels = model[["classnames"]])
  classScores <- predict(model, as.matrix(test), s = lambda, type = "response")[, , 1]
  
  if(class(classScores) == "matrix")
    classScores <- classScores[, model[["classnames"]]]
  else # Leave-one-out cross-validation likely used and glmnet doesn't have consistent return types.
    classScores <- t(classScores[model[["classnames"]]])
  
  switch(returnType, class = classPredictions, # Factor vector.
         score = classScores, # Numeric matrix.
         both = data.frame(class = classPredictions, classScores))
})

# One or more omics data sets, possibly with clinical data.
setMethod("elasticNetGLMpredictInterface", c("multnet", "MultiAssayExperiment"),
          function(model, test, targets = names(test), ...)
{
  tablesAndClasses <- .MAEtoWideTable(test, targets)
  test <- tablesAndClasses[["dataTable"]]

  elasticNetGLMpredictInterface(model, test, ...)
})