setGeneric("elasticNetGLMtrainInterface", function(measurements, ...)
{standardGeneric("elasticNetGLMtrainInterface")})

setMethod("elasticNetGLMtrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  .elasticNetGLMtrainInterface(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

# Clinical data only.
setMethod("elasticNetGLMtrainInterface", "DataFrame", function(measurements, classes, ...)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]

  .elasticNetGLMtrainInterface(measurements, splitDataset[["classes"]], ...)
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
    .elasticNetGLMtrainInterface(measurements, classes, ...)
})

.elasticNetGLMtrainInterface <- function(measurements, classes, ..., verbose = 3)
{
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting elastic net regularised GLM classifier to data.")

  glmnet::glmnet(as.matrix(measurements), classes, family = "multinomial", ...)
}

# Matrix of numeric measurements.
setGeneric("elasticNetGLMpredictInterface", function(model, test, ...)
{standardGeneric("elasticNetGLMpredictInterface")})

setMethod("elasticNetGLMpredictInterface", c("multnet", "matrix"),
          function(model, test, ...)
{
  .elasticNetGLMpredictInterface(model, DataFrame(t(test), check.names = FALSE), ...)
})

# Clinical data only.
setMethod("elasticNetGLMpredictInterface", c("multnet", "DataFrame"), function(model, test, ...)
{
  splitDataset <- .splitDataAndClasses(test, classes)
  test <- splitDataset[["measurements"]]
  
  .elasticNetGLMpredictInterface(model, test, ...)
})

# One or more omics datasets, possibly with clinical data.
setMethod("elasticNetGLMpredictInterface", c("multnet", "MultiAssayExperiment"),
          function(model, test, targets = names(test), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  test <- tablesAndClasses[["dataTable"]]
            
  .elasticNetGLMpredictInterface(model, test, ...)
})

.elasticNetGLMpredictInterface <- function(model, test, verbose = 3)
{
  if(!requireNamespace("sparsediscrim", quietly = TRUE))
    stop("The package 'sparsediscrim' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained elastic net regularised GLM classifier.")
  
  predict(model, as.matrix(test))
}