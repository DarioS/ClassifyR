setGeneric("NSCselectionInterface", function(measurements, ...)
{standardGeneric("NSCselectionInterface")})

setMethod("NSCselectionInterface", "matrix", function(measurements, classes, ...)
{
  NSCselectionInterface(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("NSCselectionInterface", "DataFrame", # Clinical data only.
          function(measurements, classes, datasetName, trained, ...,
                   selectionName = "Shrunken Centroids", verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")

  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")
  
  minError <- min(trained[["errors"]])
  threshold <- trained[["threshold"]][max(which(trained[["errors"]] == minError))]

  params <- c(list(trained), list(list(x = t(as.matrix(measurements)), y = classes, geneid = 1:ncol(measurements))), threshold)
  chosen <- as.numeric(do.call(pamr::pamr.listgenes, params)[, 1])
  
  if(is.null(S4Vectors::mcols(measurements)))
    chosen <- colnames(measurements)[chosen]
  else
    chosen <- S4Vectors::mcols(measurements)[chosen, ]
  
  if(verbose == 3)
    message("Nearest shrunken centroid feature selection completed.")
  
  SelectResult(datasetName, selectionName, list(), list(chosen)) 
})

setMethod("NSCselectionInterface", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
            
  if(ncol(measurements) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    NSCselectionInterface(measurements, classes, ...)
})