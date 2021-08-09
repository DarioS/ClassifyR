setGeneric("NSCpredictInterface", function(trained, test, ...)
standardGeneric("NSCpredictInterface"))

setMethod("NSCpredictInterface", c("pamrtrained", "matrix"), function(trained, test, ...)
{
  NSCpredictInterface(trained, DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("NSCpredictInterface", c("pamrtrained", "DataFrame"), function(trained, test, classes = NULL, ..., returnType = c("class", "score", "both"), verbose = 3)
{
  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")
  returnType <- match.arg(returnType)
  
  if(!is.null(classes)) # Remove them.
  {
    splitDataset <- .splitDataAndClasses(test, classes) # Remove classes, if present.
    test <- splitDataset[["measurements"]]
  }
 
  minError <- min(trained[["errors"]])
  threshold <- trained[["threshold"]][max(which(trained[["errors"]] == minError))]

  test <- t(as.matrix(test))   
  classPredictions <- pamr::pamr.predict(trained, test, threshold, ...)
  classScores <- pamr::pamr.predict(trained, test, threshold, type = "posterior", ...)[, levels(trained[["y"]])]
  if(!is.matrix(classScores)) # Only one sample was predicted and pamr isn't consistent with return types.
    classScores <- t(classScores)
  
  if(verbose == 3)
    message("Nearest shrunken centroid predictions made.")
  
  switch(returnType, class = classPredictions, # Factor vector.
         score = classScores, # Numeric matrix.
         both = data.frame(class = classPredictions, classScores))
})

setMethod("NSCpredictInterface", c("pamrtrained", "MultiAssayExperiment"), function(trained, test, targets = names(test), ...)
{
  test <- .MAEtoWideTable(test, targets)[["dataTable"]] # Remove any classes, if present.
  
  if(ncol(test) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    NSCpredictInterface(trained, test, ...)
})
