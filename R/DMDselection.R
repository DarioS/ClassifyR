# Mean or Median, Three Kinds of Deviations.
setGeneric("DMDselection", function(measurements, ...)
           standardGeneric("DMDselection"))

setMethod("DMDselection", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  DMDselection(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("DMDselection", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, differences = c("both", "location", "scale"),
                   trainParams, predictParams, resubstituteParams, ..., verbose = 3)
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
  
  allClassesLocationsScales <- lapply(levels(classes), function(class)
  {
    aClassMeasurements <- measurements[which(classes == class), ]
    getLocationsAndScales(aClassMeasurements, ...)
  })
  allClassesLocations <- sapply(allClassesLocationsScales, "[[", 1)
  allClassesScales <- sapply(allClassesLocationsScales, "[[", 2)
  locationsDifferences <- apply(allClassesLocations, 1, function(locations) sum(abs(c(dist(locations)))))
  scalesDifferences <- apply(allClassesScales, 1, function(scales) sum(abs(c(dist(scales)))))
  
  divergence <- 0
  if(differences %in% c("both", "location"))
    divergence <- divergence + locationsDifferences
  if(differences %in% c("both", "scale"))
    divergence <- divergence + scalesDifferences

  orderedFeatures <- order(divergence, decreasing = TRUE)
  .pickFeatures(measurements, classes, NULL, trainParams, predictParams,
                resubstituteParams, orderedFeatures, verbose)
})

# One or more omics data sets, possibly with clinical data.
setMethod("DMDselection", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]            
  DMDselection(dataTable, classes, ...)
})