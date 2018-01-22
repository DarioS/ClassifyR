setGeneric("bartlettSelection", function(measurements, ...)
{standardGeneric("bartlettSelection")})

setMethod("bartlettSelection", "matrix", # Matrix of numeric measurements.
function(measurements, classes, ...)
{
  .bartlettSelection(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("bartlettSelection", "DataFrame", # Clinical data only.
          function(measurements, classes, ...)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  .bartlettSelection(measurements, splitDataset[["classes"]], ...)
})

# One or more omics datasets, possibly with clinical data.
setMethod("bartlettSelection", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]

  if(ncol(measurements) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    .bartlettSelection(measurements, classes, ...)
})

.bartlettSelection <- function(measurements, classes,
                               datasetName, trainParams, predictParams, resubstituteParams,
                               selectionName = "Bartlett Test", verbose = 3)
{
  if(verbose == 3)
    message("Selecting features based on Bartlett statistic.")

  pValues <- apply(measurements, 2, function(featureColumn)
    bartlett.test(featureColumn, classes)[["p.value"]])
  orderedFeatures <- order(pValues)
  
  .pickFeatures(measurements, classes, datasetName, trainParams, predictParams,
                resubstituteParams, orderedFeatures, selectionName, verbose)
}