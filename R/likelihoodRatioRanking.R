setGeneric("likelihoodRatioRanking", function(measurements, ...)
           standardGeneric("likelihoodRatioRanking"))

# Matrix of numeric measurements.
setMethod("likelihoodRatioRanking", "matrix", function(measurements, classes, ...)
{
  likelihoodRatioRanking(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("likelihoodRatioRanking", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, alternative = c(location = "different", scale = "different"),
                   ..., verbose = 3)
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
  
  order(logLikelihoodRatios)
})

# One or more omics data sets, possibly with clinical data.
setMethod("likelihoodRatioRanking", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]

  if(ncol(dataTable) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
  likelihoodRatioRanking(dataTable, classes, ...)
})