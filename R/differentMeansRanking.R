setGeneric("differentMeansRanking", function(measurements, ...)
           standardGeneric("differentMeansRanking"))

# Matrix of numeric measurements.
setMethod("differentMeansRanking", "matrix", function(measurements, classes, ...)
{
  differentMeansRanking(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

# DataFrame of numeric measurements, likely created by runTests or runTest.
setMethod("differentMeansRanking", "DataFrame",
          function(measurements, classes, verbose = 3)
{
  if(!requireNamespace("genefilter", quietly = TRUE))
    stop("The package 'genefilter' could not be found. Please install it.")

  measurementsMatrix <- t(as.matrix(measurements))
  
  if(length(levels(classes)) == 2)
    statistics <- genefilter::rowttests(measurementsMatrix, classes)[, "statistic"]
  else
    statistics <- genefilter::rowFtests(measurementsMatrix, classes)[, "statistic"]
  
  if(!is.null(S4Vectors::mcols(measurements)))
    S4Vectors::mcols(measurements)[order(statistics, decreasing = TRUE)]
  else
    colnames(measurements)[order(statistics, decreasing = TRUE)]
})

# One or more omics data sets, possibly with clinical data.
setMethod("differentMeansRanking", "MultiAssayExperiment", 
          function(measurements, targets = NULL, ...)
{
  if(is.null(targets))
    stop("'targets' must be specified but was not.")
  if(length(setdiff(targets, names(measurements))))
    stop("Some values of 'targets' are not names of 'measurements' but all must be.")                            
            
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  differentMeansRanking(measurements, classes, ...)
})