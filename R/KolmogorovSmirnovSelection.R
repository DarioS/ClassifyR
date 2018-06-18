setGeneric("KolmogorovSmirnovSelection", function(measurements, ...)
           {standardGeneric("KolmogorovSmirnovSelection")})

# Matrix of numeric measurements.
setMethod("KolmogorovSmirnovSelection", "matrix", function(measurements, classes, ...)
{
  KolmogorovSmirnovSelection(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("KolmogorovSmirnovSelection", "DataFrame", # Clinical data only.
          function(measurements, classes, datasetName,
                   trainParams, predictParams, resubstituteParams, ...,
                   selectionName = "Kolmogorov-Smirnov Test", verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")

  if(verbose == 3)
    message("Selecting features by Kolmogorov Smirnov distance.")

  oneClass <- classes == levels(classes)[1]
  otherClass <- classes == levels(classes)[2]
  KSdistance <- apply(measurements, 2, function(featureColumn)
                      stats::ks.test(featureColumn[oneClass], featureColumn[otherClass], ...)[["statistic"]])

  orderedFeatures <- order(KSdistance, decreasing = TRUE)
  .pickFeatures(measurements, classes, datasetName,
                trainParams, predictParams, resubstituteParams,
                orderedFeatures, selectionName, verbose)  
})

# One or more omics data sets, possibly with clinical data.
setMethod("KolmogorovSmirnovSelection", "MultiAssayExperiment",
function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
            
  if(ncol(dataTable) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    KolmogorovSmirnovSelection(dataTable, classes, ...)
})