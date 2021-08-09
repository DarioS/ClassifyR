setGeneric("interactorDifferences", function(measurements, ...)
           standardGeneric("interactorDifferences"))

setMethod("interactorDifferences", "matrix", # Matrix of numeric measurements.
          function(measurements, ...)
{
  interactorDifferences(DataFrame(t(measurements), check.names = FALSE), ...)
})

setMethod("interactorDifferences", "DataFrame", # Possibly mixed data types.
          function(measurements, networkSets = NULL, absolute = FALSE, verbose = 3)
{
  if(is.null(networkSets))
    stop("'networkSets' is NULL but must be provided.")

  if(verbose == 3)
    message("Calculating differences between the specified interactors.")
            
  networkIDs <- names(networkSets@sets)
  allInteractions <- data.frame(do.call(rbind, networkSets@sets), networkID = rep(networkIDs, sapply(networkSets@sets, nrow)))
  keep <- allInteractions[, 1] %in% colnames(measurements) & allInteractions[, 2] %in% colnames(measurements)
  allInteractions <- allInteractions[keep, ]
  interactorTable <- as(measurements[, allInteractions[, 1]], "matrix") # Coerce to basic matrix for calculation speed.
  otherInteractorTable <- as(measurements[, allInteractions[, 2]], "matrix")
  differences <- otherInteractorTable - interactorTable
  if(absolute == TRUE)
    differences <- abs(differences)
  differences <- DataFrame(differences)
  colnames(differences) <- paste(allInteractions[, 2], '-', allInteractions[, 1])
  S4Vectors::mcols(differences) <- DataFrame(original = factor(allInteractions[, "networkID"], levels = unique(as.character(allInteractions[, "networkID"]))))
  differences
})

setMethod("interactorDifferences", "MultiAssayExperiment", # Pick one numeric table from the data set.
          function(measurements, target = NULL, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, target)
  interactorDifferences(tablesAndClasses[["dataTable"]], ...)
})