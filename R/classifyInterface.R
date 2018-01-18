setGeneric("classifyInterface", function(measurements, ...)
{standardGeneric("classifyInterface")})

setMethod("classifyInterface", "matrix", # Matrix of integer measurements.
          function(measurements, classes, test, ...)
{
  .classifyInterface(DataFrame(t(measurements[, , drop = FALSE]), check.names = FALSE),
                     classes,
                     DataFrame(t(test[, , drop = FALSE]), check.names = FALSE), ...)
})

setMethod("classifyInterface", "DataFrame", function(measurements, classes, test, ...)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  trainingMatrix <- as.matrix(splitDataset[["measurements"]])
  isInteger <- sapply(measurements, is.integer)
  measurements <- measurements[, isInteger, drop = FALSE]
  isInteger <- sapply(test, is.integer)
  testingMatrix <- as.matrix(test[, isInteger, drop = FALSE])
            
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  .classifyInterface(trainingMatrix, splitDataset[["classes"]], testingMatrix, ...)
})

setMethod("classifyInterface", "MultiAssayExperiment",
function(measurements, test, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, "integer")
  trainingMatrix <- as.matrix(tablesAndClasses[["dataTable"]])
  classes <- tablesAndClasses[["classes"]]
  testingMatrix <- .MAEtoWideTable(test, targets, "integer")
            
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  .classifyInterface(trainingMatrix, classes, testingMatrix, ...)
})

.classifyInterface <- function(measurements, classes, test, verbose, ...)
{
  if(!requireNamespace("PoiClaClu", quietly = TRUE))
    stop("The package 'PoiClaClu' could not be found. Please install it.")
  
  PoiClaClu::Classify(as.matrix(measurements), classes, as.matrix(test), ...)
}