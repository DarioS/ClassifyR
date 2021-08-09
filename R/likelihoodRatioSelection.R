setGeneric("likelihoodRatioSelection", function(measurements, ...)
           standardGeneric("likelihoodRatioSelection"))

# Matrix of numeric measurements.
setMethod("likelihoodRatioSelection", "matrix", function(measurements, classes, ...)
{
  likelihoodRatioSelection(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("likelihoodRatioSelection", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, datasetName,
                   trainParams, predictParams, resubstituteParams,
                   alternative = c(location = "different", scale = "different"),
                   ..., selectionName = "Likelihood Ratio Test (Normal)", verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")

  if(verbose == 3)
    message("Selecting features by likelihood ratio ranking.")

  allDistribution <- getLocationsAndScales(measurements, ...)
  logLikelihoodRatios <- unlist(mapply(function(featureMeasurements, scale, location)
  sum(dnorm(featureMeasurements, scale, location, log = TRUE)),
  measurements, allDistribution[[1]], allDistribution[[2]])) -
  rowSums(sapply(levels(classes), function(class)
  {
    classMeasurements <- measurements[which(classes == class), ]
    classDistribution <- getLocationsAndScales(classMeasurements, ...)
    
    unlist(mapply(function(featureMeasurements, scale, location)
    sum(dnorm(featureMeasurements, scale, location, log = TRUE)),
    classMeasurements,
    switch(alternative[["location"]], same = allDistribution[[1]], different = classDistribution[[1]]),
    switch(alternative[["scale"]], same = allDistribution[[2]], different = classDistribution[[2]])))    
  }))
  
  orderedFeatures <- order(logLikelihoodRatios)
  .pickFeatures(measurements, classes, NULL, datasetName,
                trainParams, predictParams, resubstituteParams,
                orderedFeatures, selectionName, verbose)
})

# One or more omics data sets, possibly with clinical data.
setMethod("likelihoodRatioSelection", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]

  if(ncol(dataTable) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
  likelihoodRatioSelection(dataTable, classes, ...)
})