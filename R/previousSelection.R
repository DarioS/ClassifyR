# Automated Selection of Previously Selected Features
previousSelection <- function(measurementsTrain, classesTrain, classifyResult, minimumOverlapPercent = 80,
                              .iteration, verbose = 3)
{
  if(verbose == 3)
    message("Choosing previous features.")

  previousIDs <- chosenFeatureNames(classifyResult)[[.iteration]]
  featuresIDs <- colnames(measurementsTrain)
  if(is.character(previousIDs))
  {
    safeIDs <- intersect(make.names(previousIDs), featuresIDs)
  } else { # A data frame describing the assay and variable name of the chosen feature.
    oldSafeIDs <- paste(previousIDs[, "assay"], previousIDs[, "feature"], sep = '_')
    safeIDs <- intersect(oldSafeIDs, featuresIDs)
  }
  
  commonFeatures <- intersect(safeIDs, featuresIDs)
  overlapPercent <- length(commonFeatures) / length(safeIDs) * 100
  if(overlapPercent < minimumOverlapPercent)
    signalCondition(simpleError(paste("Number of features in common between previous and current data set is lower than", minimumOverlapPercent, "percent.")))

  match(safeIDs, featuresIDs) # Return indices, not identifiers.
}
attr(previousSelection, "name") <- "previousSelection"