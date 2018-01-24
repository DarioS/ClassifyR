setGeneric("logisticRegressionTrainInterface", function(measurements, ...)
{standardGeneric("logisticRegressionTrainInterface")})

setMethod("logisticRegressionTrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  .logisticRegressionTrainInterface(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("logisticRegressionTrainInterface", "DataFrame", function(measurements, classes, ...)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  trainingMatrix <- as.matrix(splitDataset[["measurements"]])
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]

  .logisticRegressionTrainInterface(measurements, splitDataset[["classes"]], ...)
})

setMethod("logisticRegressionTrainInterface", "MultiAssayExperiment",
function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  
  if(ncol(measurements) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    .logisticRegressionTrainInterface(measurements, classes, ...)
})

.logisticRegressionTrainInterface <- function(measurements, classes, verbose = 3)
{
  if(!requireNamespace("sparsediscrim", quietly = TRUE))
    stop("The package 'sparsediscrim' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting DLDA classifier to data.")

  sparsediscrim::dlda(as.matrix(measurements), classes)
}

setGeneric("logisticRegressionPredictInterface", function(model, test, ...)
{standardGeneric("logisticRegressionPredictInterface")})

setMethod("logisticRegressionPredictInterface", c("dlda", "matrix"),
          function(model, test, ...)
{
  .logisticRegressionPredictInterface(model, DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("logisticRegressionPredictInterface", c("dlda", "DataFrame"), function(model, test, ...)
{
  splitDataset <- .splitDataAndClasses(test, classes)
  testMatrix <- splitDataset[["measurements"]]
  isNumeric <- sapply(testMatrix, is.numeric)
  testMatrix <- testMatrix[, isNumeric, drop = FALSE]
  
  .logisticRegressionPredictInterface(model, testMatrix, ...)
})

setMethod("logisticRegressionPredictInterface", c("dlda", "MultiAssayExperiment"),
          function(model, test, targets = names(test), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  test <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
            
  if(ncol(test) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    .logisticRegressionPredictInterface(model, test, classes, ...)
})

.logisticRegressionPredictInterface <- function(model, test, classes, verbose = 3)
{
  if(!requireNamespace("sparsediscrim", quietly = TRUE))
    stop("The package 'sparsediscrim' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained DLDA classifier.")
  
  predict(model, as.matrix(test))
}