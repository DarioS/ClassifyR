setGeneric("KullbackLeiblerSelection", function(measurements, ...)
           {standardGeneric("KullbackLeiblerSelection")})

# Matrix of numeric measurements.
setMethod("KullbackLeiblerSelection", "matrix", function(measurements, classes, ...)
{
  KullbackLeiblerSelection(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("KullbackLeiblerSelection", "DataFrame", # Clinical data only.
          function(measurements, classes, datasetName,
                   trainParams, predictParams, resubstituteParams, ...,
                   selectionName = "Kullback-Leibler Divergence", verbose = 3)
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

  orderedFeatures <- order(divergence, decreasing = TRUE)
  .pickFeatures(measurements, classes, datasetName,
                trainParams, predictParams, resubstituteParams,
                orderedFeatures, selectionName, verbose)  
})

# One or more omics datasets, possibly with clinical data.
setMethod("KullbackLeiblerSelection", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]

  if(ncol(dataTable) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    KullbackLeiblerSelection(dataTable, classes, ...)
})