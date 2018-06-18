setGeneric("previousSelection", function(measurements, ...)
{standardGeneric("previousSelection")})

setMethod("previousSelection", "matrix", 
          function(measurements, ...)
{
  previousSelection(DataFrame(t(measurements), check.names = FALSE), ...)
})

# Classes is passed around because most other selection functions need it, so it is sent from
# .doSelection
setMethod("previousSelection", "DataFrame", 
          function(measurements, classes, datasetName, classifyResult, minimumOverlapPercent = 80,
                               selectionName = "Previous Selection", .iteration, verbose = 3)
{
  if(verbose == 3)
    message("Choosing previous features.")
  
  if(length(.iteration) == 1)
    previousIDs <- features(classifyResult)[[.iteration]]
  else # Resample index and fold index.
    previousIDs <- features(classifyResult)[[.iteration[[1]]]][[.iteration[[2]]]]

  if(is.character(previousIDs))
  {
    commonFeatures <- intersect(previousIDs, colnames(measurements))
    overlapPercent <- length(commonFeatures) / length(previousIDs) * 100
  } else { # A data.frame describing the data set and variable name of the chosen feature.
    keepRows <- numeric()
    varInfo <- S4Vectors::mcols(measurements) # mcols stores source information about variables.
    variable <- varInfo[, "rowname"]
    variable[is.na(variable)] <- varInfo[is.na(variable), "colname"]
    for(index in 1:length(previousIDs))
    {
      if(any(previousIDs[index, "dataset"] == varInfo[, "sourceName"] & previousIDs[index, "variable"] == variable))
        keepRows <- c(keepRows, index)
    }
    commonFeatures <- previousIDs[keepRows, ]
    overlapPercent <- nrow(commonFeatures) / nrow(previousIDs) * 100
  }
  if(overlapPercent < minimumOverlapPercent)
    signalCondition(simpleError(paste("Number of features in common between previous and current data set is lower than", minimumOverlapPercent, "percent.")))
  
  SelectResult(datasetName, selectionName, list(), list(commonFeatures)) # Ranking isn't transferred across.
})

setMethod("previousSelection", "MultiAssayExperiment", 
          function(measurements, ...)
          {
            clinicalColumns <- colnames(MultiAssayExperiment::colData(clinicalColumns))
            dataTable <- wideFormat(measurements, colDataCols = clinicalColumns, check.names = FALSE, collapse = ':')
            S4Vectors::mcols(dataTable)[, "sourceName"] <- gsub("colDataCols", "clinical", S4Vectors::mcols(dataTable)[, "sourceName"])
            previousSelection(dataTable, ...)
          })