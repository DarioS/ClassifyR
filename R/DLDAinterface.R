setGeneric("DLDAtrainInterface", function(measurements, ...)
{standardGeneric("DLDAtrainInterface")})

setMethod("DLDAtrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  .DLDAtrainInterface(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("DLDAtrainInterface", "DataFrame", function(measurements, classes, ...)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  trainingMatrix <- as.matrix(splitDataset[["measurements"]])
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]

  .DLDAtrainInterface(measurements, splitDataset[["classes"]], ...)
})

setMethod("DLDAtrainInterface", "MultiAssayExperiment",
function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  
  if(ncol(measurements) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    .DLDAtrainInterface(measurements, classes, ...)
})

.DLDAtrainInterface <- function(measurements, classes, verbose = 3)
{
  if(!requireNamespace("sparsediscrim", quietly = TRUE))
    stop("The package 'sparsediscrim' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting DLDA classifier to data.")

  sparsediscrim::dlda(as.matrix(measurements), classes)
}

setGeneric("DLDApredictInterface", function(model, test, ...)
{standardGeneric("DLDApredictInterface")})

setMethod("DLDApredictInterface", c("dlda", "matrix"),
          function(model, test, ...)
{
  .DLDApredictInterface(model, DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("DLDApredictInterface", c("dlda", "DataFrame"), function(model, test, ...)
{
  splitDataset <- .splitDataAndClasses(test, classes)
  testMatrix <- splitDataset[["measurements"]]
  isNumeric <- sapply(testMatrix, is.numeric)
  testMatrix <- testMatrix[, isNumeric, drop = FALSE]
  
  .DLDApredictInterface(model, testMatrix, ...)
})

setMethod("DLDApredictInterface", c("dlda", "MultiAssayExperiment"),
          function(model, test, targets = names(test), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  test <- tablesAndClasses[["dataTable"]]
            
  if(ncol(test) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    .DLDApredictInterface(model, test, ...)
})

.DLDApredictInterface <- function(model, test, verbose = 3)
{
  if(!requireNamespace("sparsediscrim", quietly = TRUE))
    stop("The package 'sparsediscrim' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained DLDA classifier.")
  
  predict(model, as.matrix(test))
}