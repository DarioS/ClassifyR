setGeneric("leveneRanking", function(measurements, ...)
           standardGeneric("leveneRanking"))

# Matrix of numeric measurements.
setMethod("leveneRanking", "matrix", function(measurements, classes, ...)
{
  leveneRanking(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("leveneRanking", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, verbose = 3)
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
  
  if(!is.null(S4Vectors::mcols(measurements)))
    S4Vectors::mcols(measurements)[order(pValues), ]
  else
    colnames(measurements)[order(pValues)]
})

# One or more omics data sets, possibly with clinical data.
setMethod("leveneRanking", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), classes, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]

  if(ncol(dataTable) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    leveneRanking(dataTable, classes, ...)
})