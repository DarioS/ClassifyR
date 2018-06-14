setGeneric("interactorDifferences", function(measurements, ...)
           {standardGeneric("interactorDifferences")})

setMethod("interactorDifferences", "matrix", # Matrix of numeric measurements.
          function(measurements, ...)
{
  interactorDifferences(DataFrame(t(measurements), check.names = FALSE), ...)
})

setMethod("interactorDifferences", "DataFrame", # Possibly mixed data types.
          function(measurements, networkSets = NULL, verbose = 3)
{
  if(is.null(networkSets))
    stop("'networkSets' is NULL but must be provided.")

  if(verbose == 3)
    message("Calculating differences between the specified interactors.")
            
  networkIDs <- names(networkSets@sets)
  allInteractions <- do.call(rbind, networkSets@sets)
  interactorTable <- as(measurements[, allInteractions[, 1]], "matrix") # Coerce to basic matrix for calculation speed.
  otherInteractorTable <- as(measurements[, allInteractions[, 2]], "matrix")
  differences <- DataFrame(otherInteractorTable - interactorTable)
  colnames(differences) <- paste(allInteractions[, 2], '-', allInteractions[, 1])
  S4Vectors::mcols(differences) <- DataFrame(original = factor(rep(networkIDs, sapply(networkSets@sets, nrow)), levels = networkIDs))
  differences
})

setMethod("interactorDifferences", "MultiAssayExperiment", # Pick one numeric table from the data set.
          function(measurements, target = NULL, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, target)
  interactorDifferences(tablesAndClasses[["dataTable"]], ...)
})