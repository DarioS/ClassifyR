setGeneric("leveneSelection", function(measurements, ...)
           {standardGeneric("leveneSelection")})

# Matrix of numeric measurements.
setMethod("leveneSelection", "matrix", function(measurements, classes, ...)
{
  leveneSelection(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("leveneSelection", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, datasetName,
                   trainParams, predictParams, resubstituteParams,
                   selectionName = "Levene Test", verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  
  if(!requireNamespace("car", quietly = TRUE))
    stop("The package 'car' could not be found. Please install it.")            
  if(verbose == 3)
    message("Calculating Levene statistic.")

  pValues <- apply(measurements, 2, function(featureColumn)
             car::leveneTest(featureColumn, classes)[["Pr(>F)"]][1])
  orderedFeatures <- order(pValues)
 
  .pickFeatures(measurements, classes, NULL, datasetName,
                trainParams, predictParams, resubstituteParams,
                orderedFeatures, selectionName, verbose)  
})

# One or more omics data sets, possibly with clinical data.
setMethod("leveneSelection", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]

  if(ncol(dataTable) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    leveneSelection(dataTable, classes, ...)
})