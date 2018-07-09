setGeneric("kNNinterface", function(measurements, ...) {standardGeneric("kNNinterface")})

setMethod("kNNinterface", "matrix",
          function(measurements, classes, test, ...)
{
  kNNinterface(DataFrame(t(measurements), check.names = FALSE),
               classes,
               DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("kNNinterface", "DataFrame", function(measurements, classes, test, ..., verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  trainingMatrix <- as.matrix(splitDataset[["measurements"]])
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  isNumeric <- sapply(test, is.numeric)
  test <- test[, isNumeric, drop = FALSE]  
  
  if(!requireNamespace("class", quietly = TRUE))
    stop("The package 'class' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting k Nearest Neighbours classifier to data and predicting classes.")
  
  class::knn(as.matrix(measurements), as.matrix(test), classes, ...)
})

setMethod("kNNinterface", "MultiAssayExperiment",
function(measurements, test, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  trainingTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  testingTable <- .MAEtoWideTable(test, targets)
            
  .checkVariablesAndSame(trainingTable, testingTable)
  kNNinterface(trainingTable, classes, testingTable, ...)
})