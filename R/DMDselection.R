# Mean or Median, Three Kinds of Deviations.
setGeneric("DMDselection", function(measurements, ...)
           {standardGeneric("DMDselection")})

setMethod("DMDselection", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  .DMDselection(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("DMDselection", "DataFrame", # Clinical data only.
          function(measurements, classes, ...)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- apply(measurements, 2, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  .DMDselection(measurements, splitDataset[["classes"]], ...)
})

# One or more omics datasets, possibly with clinical data.
setMethod("DMDselection", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
{
  tableAndClasses <- .MAEtoWideTable(measurements, targets)
  dataTable <- tableAndClasses[["dataTable"]]
  classes <- tableAndClasses[["classes"]]            
  isNumeric <- apply(dataTable, 2, is.numeric)
  dataTable <- dataTable[, isNumeric, drop = FALSE]
  if(ncol(dataTable) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    .DMDselection(dataTable, classes, ...)
})

.DMDselection <- function(measurements, classes,
                          datasetName, trainParams, predictParams, resubstituteParams, ...,
                          selectionName = "Differences of Medians and Deviations", verbose = 3)
{
  if(verbose == 3)
    message("Selecting features by DMD.")
  
  oneClassTraining <- which(classes == levels(classes)[1])
  otherClassTraining <- which(classes == levels(classes)[2])
  oneClassMeasurements <- measurements[oneClassTraining, ]
  otherClassMeasurements <- measurements[otherClassTraining, ]
  oneClassDistribution <- getLocationsAndScales(oneClassMeasurements, ...)
  otherClassDistribution <- getLocationsAndScales(otherClassMeasurements, ...)
  locationDifference <- oneClassDistribution[[1]] - otherClassDistribution[[1]]
  divergence <- abs(locationDifference) + abs(oneClassDistribution[[2]] - otherClassDistribution[[2]])
  orderedFeatures <- order(divergence, decreasing = TRUE)

  .pickFeatures(measurements, classes, datasetName, trainParams, predictParams,
                resubstituteParams, orderedFeatures, selectionName, verbose)
}