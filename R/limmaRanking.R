setGeneric("limmaRanking", function(measurements, ...)
           standardGeneric("limmaRanking"))

# Matrix of numeric measurements.
setMethod("limmaRanking", "matrix", function(measurements, classes, ...)
{
  limmaRanking(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

# DataFrame of numeric measurements, likely created by runTests or runTest.
setMethod("limmaRanking", "DataFrame",
          function(measurements, classes, ..., verbose = 3)
{
  if(!requireNamespace("limma", quietly = TRUE))
    stop("The package 'limma' could not be found. Please install it.")

  fitParams <- list(t(as.matrix(measurements)), model.matrix(~ classes))
  if(!missing(...))
    fitParams <- append(fitParams, ...)
  linearModel <- do.call(limma::lmFit, fitParams)
  linearModel <- limma::eBayes(linearModel)
  match(rownames(limma::topTable(linearModel, 2, number = Inf, sort.by = "p")), colnames(measurements))
})

# One or more omics data sets, possibly with clinical data.
setMethod("limmaRanking", "MultiAssayExperiment", 
          function(measurements, targets = NULL, ...)
{
  if(is.null(targets))
    stop("'targets' must be specified but was not.")
  if(length(setdiff(targets, names(measurements))))
    stop("Some values of 'targets' are not names of 'measurements' but all must be.")                            
            
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  limmaRanking(measurements, classes, ...)
})