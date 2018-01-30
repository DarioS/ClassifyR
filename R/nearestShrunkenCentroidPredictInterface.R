setGeneric("nearestShrunkenCentroidPredictInterface", function(trained, test, ...)
{standardGeneric("nearestShrunkenCentroidPredictInterface")})

setMethod("nearestShrunkenCentroidPredictInterface", c("pamrtrained", "matrix"), function(trained, test, ...)
{
  nearestShrunkenCentroidPredictInterface(trained, DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("nearestShrunkenCentroidPredictInterface", c("pamrtrained", "DataFrame"), function(trained, test, ..., verbose = 3)
{
  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")
 
  minError <- min(trained[["errors"]])
  threshold <- trained[["threshold"]][max(which(trained[["errors"]] == minError))]
  
  test <- t(as.matrix(test))   
  predictions <- pamr::pamr.predict(trained, test, threshold, ...)
  
  if(verbose == 3)
    message("Nearest shrunken centroid predictions made.")
  predictions
})

setMethod("nearestShrunkenCentroidPredictInterface", c("pamrtrained", "MultiAssayExperiment"), function(trained, test, targets = names(test), ...)
{
  test <- .MAEtoWideTable(test, targets)[["dataTable"]]
  
  if(ncol(test) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    nearestShrunkenCentroidPredictInterface(trained, test, classes, ...)
})