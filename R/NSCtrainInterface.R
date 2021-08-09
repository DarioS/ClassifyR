setGeneric("NSCtrainInterface", function(measurements, ...)
standardGeneric("NSCtrainInterface"))

setMethod("NSCtrainInterface", "matrix", function(measurements, classes, ...)
{
  NSCtrainInterface(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("NSCtrainInterface", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, ..., verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")

  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")

  trainedModel <- pamr::pamr.train(list(x = t(as.matrix(measurements)), y = classes), ...)
  
  if(verbose == 3)
    message("Nearest shrunken centroid training completed.")
  
  trainedModel  
})

setMethod("NSCtrainInterface", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
{ 
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  
  if(ncol(measurements) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    NSCtrainInterface(measurements, classes, ...)
})