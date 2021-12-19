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
    pValues <- genefilter::rowttests(measurementsMatrix, classes)[, "p.value"]
  else
    pValues <- genefilter::rowFtests(measurementsMatrix, classes)[, "p.value"]
  
  if(!is.null(S4Vectors::mcols(measurements)))
    S4Vectors::mcols(measurements)[order(pValues), ]
  else
    colnames(measurements)[order(pValues)]
})

# One or more omics data sets, possibly with clinical data.
setMethod("differentMeansRanking", "MultiAssayExperiment", 
          function(measurements, targets = NULL, classes, ...)
{
  if(is.null(targets))
    stop("'targets' must be specified but was not.")
  if(length(setdiff(targets, names(measurements))))
    stop("Some values of 'targets' are not names of 'measurements' but all must be.")                            
            
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  differentMeansRanking(measurements, classes, ...)
})