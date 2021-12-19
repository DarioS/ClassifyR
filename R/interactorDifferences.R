setGeneric("interactorDifferences", function(measurements, ...)
           standardGeneric("interactorDifferences"))

setMethod("interactorDifferences", "matrix", # Matrix of numeric measurements.
          function(measurements, ...)
{
  interactorDifferences(DataFrame(t(measurements), check.names = FALSE), ...)
})

setMethod("interactorDifferences", "DataFrame", # Possibly mixed data types.
          function(measurements, featurePairs = NULL, absolute = FALSE, verbose = 3)
{
  if(is.null(featurePairs))
    stop("'featurePairs' is NULL but must be provided.")

  if(verbose == 3)
    message("Calculating differences between the specified interactors.")
            
  keep <- S4Vectors::first(featurePairs) %in% colnames(measurements) & S4Vectors::second(featurePairs) %in% colnames(measurements)
  featurePairs <- featurePairs[keep]
  interactorTable <- as(measurements[, S4Vectors::first(featurePairs)], "matrix") # Coerce to basic matrix for calculation speed.
  otherInteractorTable <- as(measurements[, S4Vectors::second(featurePairs)], "matrix")
  differences <- otherInteractorTable - interactorTable
  if(absolute == TRUE)
    differences <- abs(differences)
  differences <- DataFrame(differences)
  colnames(differences) <- paste(S4Vectors::second(featurePairs), '-', S4Vectors::first(featurePairs))
  differences
})

setMethod("interactorDifferences", "MultiAssayExperiment", # Pick one numeric table from the data set.
          function(measurements, target = NULL, classes, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, target, classes)
  interactorDifferences(tablesAndClasses[["dataTable"]], ...)
})