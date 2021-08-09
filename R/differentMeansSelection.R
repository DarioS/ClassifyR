setGeneric("differentMeansSelection", function(measurements, ...)
           standardGeneric("differentMeansSelection"))

# Matrix of numeric measurements.
setMethod("differentMeansSelection", "matrix", function(measurements, classes, ...)
{
  differentMeansSelection(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

# DataFrame of numeric measurements, likely created by runTests or runTest.
setMethod("differentMeansSelection", "DataFrame",
          function(measurements, classes, datasetName,
                   trainParams, predictParams, resubstituteParams,
                   selectionName = "Difference in Means", verbose = 3)
{
  if(!requireNamespace("genefilter", quietly = TRUE))
    stop("The package 'genefilter' could not be found. Please install it.")

  measurementsMatrix <- t(as.matrix(measurements))
  
  if(length(levels(classes)) == 2)
    statistics <- genefilter::rowttests(measurementsMatrix, classes)[, "statistic"]
  else
    statistics <- genefilter::rowFtests(measurementsMatrix, classes)[, "statistic"]
  orderedFeatures <- order(statistics, decreasing = TRUE)

  .pickFeatures(measurements, classes, NULL,
                datasetName, trainParams, predictParams, resubstituteParams,
                orderedFeatures, selectionName, verbose)  
})

# One or more omics data sets, possibly with clinical data.
setMethod("differentMeansSelection", "MultiAssayExperiment", 
          function(measurements, targets = NULL, ...)
{
  if(is.null(targets))
    stop("'targets' must be specified but was not.")
  if(length(setdiff(targets, names(measurements))))
    stop("Some values of 'targets' are not names of 'measurements' but all must be.")                            
            
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  differentMeansSelection(measurements, classes, ...)
})