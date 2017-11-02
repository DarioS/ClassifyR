setGeneric("nearestShrunkenCentroidPredictInterface", function(trained, test, ...)
{standardGeneric("nearestShrunkenCentroidPredictInterface")})

setMethod("nearestShrunkenCentroidPredictInterface", c("pamrtrained", "matrix"), function(trained, test, ...)
{
  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")
  
  nearestShrunkenCentroidPredictInterface(trained, ExpressionSet(test), ...)
})

setMethod("nearestShrunkenCentroidPredictInterface", c("pamrtrained", "ExpressionSet"), function(trained, test, ..., verbose = 3)
{
  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")
  test <- exprs(test)
  
  minError <- min(trained[["errors"]])
  threshold <- trained[["threshold"]][max(which(trained[["errors"]] == minError))]
  
  predictions <- pamr::pamr.predict(trained, test, threshold, ...)
  
  if(verbose == 3)
    message("Nearest shrunken centroid predictions made.")
  predictions
})