# Mean or Median, Three Kinds of Deviations.
setGeneric("DMDselection", function(measurements, ...)
           {standardGeneric("DMDselection")})

setMethod("DMDselection", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  DMDselection(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("DMDselection", "DataFrame", # Clinical data only.
          function(measurements, classes, datasetName, differences = c("both", "location", "scale"),
                   trainParams, predictParams, resubstituteParams, ...,
                   selectionName = "Differences of Medians and Deviations", verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  
  if(verbose == 3)
    message("Selecting features by DMD.")
  differences <- match.arg(differences)
  
  oneClassTraining <- which(classes == levels(classes)[1])
  otherClassTraining <- which(classes == levels(classes)[2])
  oneClassMeasurements <- measurements[oneClassTraining, ]
  otherClassMeasurements <- measurements[otherClassTraining, ]
  oneClassDistribution <- getLocationsAndScales(oneClassMeasurements, ...)
  otherClassDistribution <- getLocationsAndScales(otherClassMeasurements, ...)
  locationDifference <- abs(oneClassDistribution[[1]] - otherClassDistribution[[1]])
  scaleDifference <- abs(oneClassDistribution[[2]] - otherClassDistribution[[2]])
  
  divergence <- 0
  if(differences %in% c("both", "location"))
    divergence <- divergence + locationDifference
  if(differences %in% c("both", "scale"))
    divergence <- divergence + scaleDifference

  orderedFeatures <- order(divergence, decreasing = TRUE)
  .pickFeatures(measurements, classes, datasetName, trainParams, predictParams,
                resubstituteParams, orderedFeatures, selectionName, verbose)  
})

# One or more omics datasets, possibly with clinical data.
setMethod("DMDselection", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]            
  DMDselection(dataTable, classes, ...)
})