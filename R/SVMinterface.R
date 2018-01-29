setGeneric("SVMtrainInterface", function(measurements, ...)
{standardGeneric("SVMtrainInterface")})

setMethod("SVMtrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  .SVMtrainInterface(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("SVMtrainInterface", "DataFrame", function(measurements, classes, ...)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  trainingMatrix <- as.matrix(splitDataset[["measurements"]])
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]

  .SVMtrainInterface(measurements, splitDataset[["classes"]], ...)
})

setMethod("SVMtrainInterface", "MultiAssayExperiment",
function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  
  if(ncol(measurements) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    .SVMtrainInterface(measurements, classes, ...)
})

.SVMtrainInterface <- function(measurements, classes, ..., verbose = 3)
{
  if(!requireNamespace("e1071", quietly = TRUE))
    stop("The package 'e1071' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting SVM classifier to data.")

  svm(measurements, classes, ...)
}

setGeneric("SVMpredictInterface", function(model, test, ...)
{standardGeneric("SVMpredictInterface")})

setMethod("SVMpredictInterface", c("svm", "matrix"),
          function(model, test, ...)
{
  .SVMpredictInterface(model, DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("SVMpredictInterface", c("svm", "DataFrame"), function(model, test, ...)
{
  splitDataset <- .splitDataAndClasses(test, classes)
  testMatrix <- splitDataset[["measurements"]]
  isNumeric <- sapply(testMatrix, is.numeric)
  testMatrix <- testMatrix[, isNumeric, drop = FALSE]
  
  .SVMpredictInterface(model, testMatrix, ...)
})

setMethod("SVMpredictInterface", c("svm", "MultiAssayExperiment"),
          function(model, test, targets = names(test), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  test <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
            
  if(ncol(test) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    .SVMpredictInterface(model, test, classes, ...)
})

.SVMpredictInterface <- function(model, test, classes, verbose = 3)
{
  if(!requireNamespace("sparsediscrim", quietly = TRUE))
    stop("The package 'sparsediscrim' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained SVM classifier.")
  
  predict(model, test)
}