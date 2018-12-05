setGeneric("elasticNetFeatures", function(model, ...)
           {standardGeneric("elasticNetFeatures")})

setMethod("elasticNetFeatures", "multnet",
          function(model)
{
  inputFeatures <- attr(model, "features")
  browser()          
  # Floating point numbers test for equality.
  whichCoefficientColumn <- which(abs(model[["lambda"]] - attr(model, "tune")[["lambda"]]) < 0.00001)[1]
  coefficientsUsed <- sapply(model[["beta"]], function(classCoefficients) classCoefficients[, whichCoefficientColumn])
  featureScores <- rowSums(abs(coefficientsUsed))
  rankedFeatures <- inputFeatures[order(featureScores, decreasing = TRUE)]
  selectedFeatures <- inputFeatures[featureScores != 0]
  
  # Colon is a reserved symbol for separating data name and feature name, which is necessary for identifiability of MultiAssayExperiment features. It is not permitted in feature names.
  if(grepl(':', selectedFeatures[1]) == TRUE) # Convert to data.frame.
  {
    selectedFeatures <- do.call(rbind, strsplit(selectedFeatures, ':'))
    rankedFeatures <- do.call(rbind, strsplit(rankedFeatures, ':'))
    colnames(selectedFeatures) <- colnames(rankedFeatures) <- c("dataset", "feature")
  }
  list(unique(rankedFeatures), selectedFeatures)
})