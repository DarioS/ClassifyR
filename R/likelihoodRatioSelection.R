setGeneric("likelihoodRatioSelection", function(measurements, ...)
           {standardGeneric("likelihoodRatioSelection")})

# Matrix of numeric measurements.
setMethod("likelihoodRatioSelection", "matrix", function(measurements, classes, ...)
{
  likelihoodRatioSelection(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("likelihoodRatioSelection", "DataFrame", # Clinical data only.
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
  
  oneClass <- classes == levels(classes)[1]
  otherClass <- classes == levels(classes)[2]
  oneClassMeasurements <- measurements[oneClass, ]
  otherClassMeasurements <- measurements[otherClass, ]
  oneClassDistribution <- getLocationsAndScales(oneClassMeasurements, ...)
  otherClassDistribution <- getLocationsAndScales(otherClassMeasurements, ...)
  allDistribution <- getLocationsAndScales(measurements, ...)

  logLikelihoodRatios <- -2 * (unlist(mapply(function(featureMeasurements, scale, location)
  sum(dnorm(featureMeasurements, scale, location, log = TRUE)),
  measurements, allDistribution[[1]], allDistribution[[2]])) -
  unlist(mapply(function(featureMeasurements, scale, location)
  sum(dnorm(featureMeasurements, scale, location, log = TRUE)),
  oneClassMeasurements,
  switch(alternative[["location"]], same = allDistribution[[1]], different = oneClassDistribution[[1]]),
  switch(alternative[["scale"]], same = allDistribution[[2]], different = oneClassDistribution[[2]]))) -
  unlist(mapply(function(featureMeasurements, scale, location)
  sum(dnorm(featureMeasurements, scale, location, log = TRUE)),
  otherClassMeasurements,
  switch(alternative[["location"]], same = allDistribution[[1]], different = otherClassDistribution[[1]]),
  switch(alternative[["scale"]], same = allDistribution[[2]], different = otherClassDistribution[[2]]))))
  orderedFeatures <- order(logLikelihoodRatios, decreasing = TRUE)
  
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