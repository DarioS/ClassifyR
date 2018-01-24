setGeneric("nearestShrunkenCentroidSelectionInterface", function(measurements, ...)
{standardGeneric("nearestShrunkenCentroidSelectionInterface")})

setMethod("nearestShrunkenCentroidSelectionInterface", "matrix", function(measurements, classes, ...)
{
  .nearestShrunkenCentroidSelectionInterface(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("nearestShrunkenCentroidSelectionInterface", "DataFrame", # Clinical data only.
          function(measurements, classes, ...)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  .nearestShrunkenCentroidSelectionInterface(measurements, splitDataset[["classes"]], ...)
})

setMethod("nearestShrunkenCentroidSelectionInterface", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
            
  if(ncol(measurements) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    .nearestShrunkenCentroidSelectionInterface(measurements, classes, ...)
})

.nearestShrunkenCentroidSelectionInterface <- function(measurements, classes, datasetName, trained, ..., selectionName = "Shrunken Centroids",
                                                      verbose = 3)
{
  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")
  
  minError <- min(trained[["errors"]])
  threshold <- trained[["threshold"]][max(which(trained[["errors"]] == minError))]
  
  params <- list(...)
  params <- params[!names(params) %in% c("trainParams", "predictParams")]
  params <- c(list(trained), list(list(x = t(as.matrix(measurements)), y = classes, geneid = 1:ncol(measurements))), threshold, params)

  chosen <- as.numeric(do.call(pamr::pamr.listgenes, params)[, 1])
  
  if(is.null(mcols(measurements)))
    chosen <- colnames(measurements)[chosen]
  else
    chosen <- mcols(measurements)[chosen, ]
  
  if(verbose == 3)
    message("Nearest shrunken centroid feature selection completed.")
  
  SelectResult(datasetName, selectionName, list(), list(chosen))
}