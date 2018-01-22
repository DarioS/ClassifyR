setGeneric("previousSelection", function(measurements, ...)
{standardGeneric("previousSelection")})

setMethod("previousSelection", "matrix", 
          function(measurements, ...)
          {
            .previousSelection(DataFrame(t(measurements), check.names = FALSE), ...)
          })

setMethod("previousSelection", "DataFrame", 
          function(measurements, ...)
          {
            .previousSelection(measurements, ...)
          })

setMethod("previousSelection", "MultiAssayExperiment", 
          function(measurements, ...)
          {
            clinicalColumns <- colnames(colData(clinicalColumns))
            dataTable <- wideFormat(measurements, colDataCols = clinicalColumns, check.names = FALSE)
            .previousSelection(dataTable, ...)
          })

.previousSelection <- function(measurements, datasetName, classifyResult, minimumOverlapPercent = 80,
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
  } else { # A data.frame describing the dataset and variable name of the chosen feature.
    keepRows <- numeric()
    varInfo <- mcols(measurements)
    for(index in 1:length(previousIDs)) # mcols stores source information about variables.
    {
      if(any(previousIDs[1, "dataset"] %in% varInfo[, "dataset"] & previousIDs[1, "variable"] %in% varInfo[, "variable"]))
        keepRows <- c(keepRows, index)
    }
    commonFeatures <- previousIDs[keepRows, ]
    overlapPercent <- nrow(commonFeatures) / nrow(previousIDs) * 100
  }
  if(overlapPercent < minimumOverlapPercent)
    signalCondition(simpleError(paste("Number of features in common between previous and current dataset is lower than", minimumOverlapPercent, "percent.")))
  
  SelectResult(datasetName, selectionName, list(), list(commonFeatures)) # Ranking isn't transferred across.
}