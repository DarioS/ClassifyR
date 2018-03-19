setGeneric("elasticNetGLMtrainInterface", function(measurements, ...)
{standardGeneric("elasticNetGLMtrainInterface")})

setMethod("elasticNetGLMtrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  elasticNetGLMtrainInterface(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

# Clinical data only.
setMethod("elasticNetGLMtrainInterface", "DataFrame", function(measurements, classes, lambda = NULL, ..., verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]

  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting elastic net regularised GLM classifier to data.")

  fitted <- glmnet::glmnet(as.matrix(measurements), splitDataset[["classes"]], family = "multinomial", ...)
  
  if(is.null(lambda)) # fitted has numerous models for automatically chosen lambda values.
  { # Pick one lambda based on resubstitution performance.
    bestLambda <- fitted[["lambda"]][which.max(fitted[["dev.ratio"]])[1]]
    attr(fitted, "tune") <- list(lambda = bestLambda)
  }
  
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
setMethod("elasticNetGLMpredictInterface", c("multnet", "DataFrame"), function(model, test, classes = NULL, lambda, ..., verbose = 3)
{ # ... just consumes emitted tuning variables from .doTrain which are unused.
  if(!is.null(classes))
  {
    splitDataset <- .splitDataAndClasses(test, classes)  # Remove any classes, if present.
    test <- splitDataset[["measurements"]]
  }
  
  if(!requireNamespace("sparsediscrim", quietly = TRUE))
    stop("The package 'sparsediscrim' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained elastic net regularised GLM classifier.")

  if(missing(lambda)) # Tuning parameters are not passed to prediction functions.
    lambda <- attr(model, "tune")[["lambda"]] # Sneak it in as an attribute on the model.
  factor(predict(model, as.matrix(test), s = lambda, type = "class"))
})

# One or more omics datasets, possibly with clinical data.
setMethod("elasticNetGLMpredictInterface", c("multnet", "MultiAssayExperiment"),
          function(model, test, targets = names(test), ...)
{
  tablesAndClasses <- .MAEtoWideTable(test, targets)
  test <- tablesAndClasses[["dataTable"]]

  elasticNetGLMpredictInterface(model, test, ...)
})