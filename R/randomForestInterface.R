setGeneric("randomForestInterface", function(measurements, ...)
{standardGeneric("randomForestInterface")})

setMethod("randomForestInterface", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, test, ...)
{
  randomForestInterface(DataFrame(t(measurements), check.names = FALSE),
                        classes,
                        DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("randomForestInterface", "DataFrame", function(measurements, classes, test, ..., verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  test <- .splitDataAndClasses(test, classes)[["measurements"]]

  if(!requireNamespace("randomForest", quietly = TRUE))
    stop("The package 'randomForestInterface' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting random forest classifier to training data and making predictions on test
            data.")

  # Convert to base data.frame as randomForest doesn't understand DataFrame.
  randomForest::randomForest(as.data.frame(splitDataset[["measurements"]]), splitDataset[["classes"]], as.data.frame(test), ...)
})

setMethod("randomForestInterface", "MultiAssayExperiment",
function(measurements, targets = names(measurements), test, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, restrict = NULL)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  test <- .MAEtoWideTable(measurements, targets, restrict = NULL)[["measurements"]]
  
  randomForestInterface(measurements, classes, test, ...)
})