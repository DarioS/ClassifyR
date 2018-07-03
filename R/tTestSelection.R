setGeneric("tTestSelection", function(measurements, ...)
           {standardGeneric("tTestSelection")})

# Matrix of numeric measurements.
setMethod("tTestSelection", "matrix", function(measurements, classes, ...)
{
  tTestSelection(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

# DataFrame of numeric measurements, likely created by runTests or runTest.
setMethod("tTestSelection", "DataFrame",
          function(measurements, classes, datasetName,
                   trainParams, predictParams, resubstituteParams,
                   selectionName = "t-test", verbose = 3)
{
  if(!requireNamespace("genefilter", quietly = TRUE))
    stop("The package 'genefilter' could not be found. Please install it.")

  measurementsMatrix <- t(as.matrix(measurements))
  tStats <- genefilter::rowttests(measurementsMatrix, classes)[, "statistic"]
  orderedFeatures <- order(tStats, decreasing = TRUE)

  .pickFeatures(measurements, classes, NULL,
                datasetName, trainParams, predictParams, resubstituteParams,
                orderedFeatures, selectionName, verbose)  
})

# One or more omics data sets, possibly with clinical data.
setMethod("tTestSelection", "MultiAssayExperiment", 
          function(measurements, targets = NULL, ...)
{
  if(is.null(targets))
    stop("'targets' must be specified but was not.")
  if(length(setdiff(targets, names(measurements))))
    stop("Some values of 'targets' are not names of 'measurements' but all must be.")                            
            
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  tTestSelection(measurements, classes, ...)
})