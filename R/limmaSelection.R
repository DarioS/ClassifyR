setGeneric("limmaSelection", function(measurements, ...)
           {standardGeneric("limmaSelection")})

# Matrix of numeric measurements.
setMethod("limmaSelection", "matrix", function(measurements, classes, ...)
{
  .limmaSelection(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

# DataFrame of numeric measurements, likely created by runTests or runTest.
setMethod("limmaSelection", "DataFrame", function(measurements, classes, ...)
{
  .limmaSelection(measurements, classes, ...)
})

# One or more omics datasets, possibly with clinical data.
setMethod("limmaSelection", "MultiAssayExperiment", 
          function(measurements, targets = NULL, ...)
{
  if(is.null(targets))
    stop("'targets' must be specified but was not.")
  if(length(setdiff(targets, names(measurements))))
    stop("Some values of 'targets' are not names of 'measurements' but all must be.")                            
            
  if(verbose == 3)
    message("Doing feature selection.")
            
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  .limmaSelection(measurements, classes, ...)
})

.limmaSelection <- function(measurements, classes, datasetName,
                            trainParams, predictParams, resubstituteParams, ...,
                            selectionName = "Moderated t-test", verbose = 3)
{
  if(!requireNamespace("limma", quietly = TRUE))
    stop("The package 'limma' could not be found. Please install it.")

  fitParams <- list(t(as.matrix(measurements)), model.matrix(~ classes))
  if(!missing(...))
    fitParams <- append(fitParams, ...)
  linearModel <- do.call(limma::lmFit, fitParams)
  linearModel <- limma::eBayes(linearModel)
  orderedFeatures <- match(rownames(limma::topTable(linearModel, 2, number = Inf, sort.by = "p")),
                           colnames(measurements))

  .pickFeatures(measurements, classes,
                datasetName, trainParams, predictParams, resubstituteParams,
                orderedFeatures, selectionName, verbose)  
}