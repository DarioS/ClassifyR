setGeneric("leveneSelection", function(measurements, ...)
           {standardGeneric("leveneSelection")})

# Matrix of numeric measurements.
setMethod("leveneSelection", "matrix", function(measurements, classes, ...)
{
  .leveneSelection(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("leveneSelection", "DataFrame", # Clinical data only.
          function(measurements, classes, ...)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  .leveneSelection(measurements, splitDataset[["classes"]], ...)
})

# One or more omics datasets, possibly with clinical data.
setMethod("leveneSelection", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]

  if(ncol(dataTable) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    .leveneSelection(dataTable, classes, ...)
})

.leveneSelection <- function(measurements, classes, datasetName,
                             trainParams, predictParams, resubstituteParams,
                             selectionName = "Levene Test", verbose = 3)
{
  if(!requireNamespace("car", quietly = TRUE))
    stop("The package 'car' could not be found. Please install it.")            
  if(verbose == 3)
    message("Calculating Levene statistic.")

  pValues <- apply(measurements, 2, function(featureColumn)
             car::leveneTest(featureColumn, classes)[["Pr(>F)"]][1])
  orderedFeatures <- order(pValues)
 
  .pickFeatures(measurements, classes, datasetName,
                trainParams, predictParams, resubstituteParams,
                orderedFeatures, selectionName, verbose)
}
