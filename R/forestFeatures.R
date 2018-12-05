setGeneric("forestFeatures", function(forest, ...)
           {standardGeneric("forestFeatures")})

setMethod("forestFeatures", "randomForest",
          function(forest)
{
  inputFeatures <- rownames(randomForest::importance(forest))
  rankedFeatures <- inputFeatures[order(randomForest::importance(forest), decreasing = TRUE)]
  selectedFeatures <- inputFeatures[randomForest::varUsed(forest) > 0]
  selectedFeatures <- selectedFeatures[na.omit(match(rankedFeatures, selectedFeatures))]
  
  # Colon is a reserved symbol for separating data name and feature name, which is necessary for identifiability of MultiAssayExperiment features. It is not permitted in feature names.
  if(grepl(':', selectedFeatures[1]) == TRUE) # Convert to data.frame.
  {
    selectedFeatures <- do.call(rbind, strsplit(selectedFeatures, ':'))
    rankedFeatures <- do.call(rbind, strsplit(rankedFeatures, ':'))
    colnames(selectedFeatures) <- colnames(rankedFeatures) <- c("dataset", "feature")
  }
  list(rankedFeatures, selectedFeatures)
})