setGeneric("KullbackLeiblerRanking", function(measurements, ...)
           standardGeneric("KullbackLeiblerRanking"))

# Matrix of numeric measurements.
setMethod("KullbackLeiblerRanking", "matrix", function(measurements, classes, ...)
{
  KullbackLeiblerRanking(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("KullbackLeiblerRanking", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, ..., verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  
  if(verbose == 3)
    message("Selecting features by Kullback-Leibler divergence.")

  oneClassMeasurements <- measurements[classes == levels(classes)[1], ]
  otherClassMeasurements <- measurements[classes == levels(classes)[2], ]
  oneClassDistribution <- getLocationsAndScales(oneClassMeasurements, ...)
  otherClassDistribution <- getLocationsAndScales(otherClassMeasurements, ...)
  locationDifference <- oneClassDistribution[[1]] - otherClassDistribution[[1]]
  divergence <- 1/2 * (locationDifference^2 / ((oneClassDistribution[[2]])^2) +
                         locationDifference^2 / ((otherClassDistribution[[2]])^2) +
                         ((oneClassDistribution[[2]])^2) / ((otherClassDistribution[[2]])^2) +
                         ((otherClassDistribution[[2]])^2) / ((oneClassDistribution[[2]])^2))

  if(!is.null(S4Vectors::mcols(measurements)))
    S4Vectors::mcols(measurements)[order(divergence, decreasing = TRUE), ]
  else
    colnames(measurements)[order(divergence, decreasing = TRUE)]
  
})

# One or more omics data sets, possibly with clinical data.
setMethod("KullbackLeiblerRanking", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), classes, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]

  if(ncol(dataTable) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    KullbackLeiblerRanking(dataTable, classes, ...)
})