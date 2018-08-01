setGeneric("kTSPclassifier", function(measurements, ...)
           {standardGeneric("kTSPclassifier")})

setMethod("kTSPclassifier", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, test, featurePairs = NULL, ...)
{
  if(is.null(featurePairs))
    stop("No feature pairs provided but some must be.")
  if(!"Pairs" %in% class(featurePairs))
    stop("'featurePairs' must be of type Pairs.")            
            
  kTSPclassifier(DataFrame(t(measurements[, , drop = FALSE]), check.names = FALSE),
                 classes,
                 DataFrame(t(test[, , drop = FALSE]), check.names = FALSE), featurePairs, ...)
})

setMethod("kTSPclassifier", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, test, featurePairs = NULL,
                   weighted = c("unweighted", "weighted", "both"), minDifference = 0,
                   returnType = c("class", "score", "both"), verbose = 3)
{
  if(is.null(featurePairs))
    stop("No feature pairs provided but some must be.")
  if(!"Pairs" %in% class(featurePairs))
    stop("'featurePairs' must be of type Pairs.")            
            
  splitDataset <- .splitDataAndClasses(measurements, classes)
  trainingMatrix <- splitDataset[["measurements"]]
  isNumeric <- sapply(trainingMatrix, is.numeric)
  trainingMatrix <- as.matrix(trainingMatrix[, isNumeric, drop = FALSE])
  isNumeric <- sapply(test, is.numeric)
  testingMatrix <- as.matrix(test[, isNumeric, drop = FALSE])
  
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  
  weighted <- match.arg(weighted)
  returnType <- match.arg(returnType)
  
  classesSizes <- sapply(levels(classes), function(class) sum(classes == class))
  largerClass <- names(classesSizes)[which.max(classesSizes)[1]]
  secondClass <- classes == levels(classes)[2]
  
  if(verbose == 3)
    message("Determining inequalities of feature pairs.")

  # Order pairs so that first < second is the rule for predicting the second class, based on factor levels.
  # Effectively the classifier training.
  featurePairs <- do.call(c, lapply(featurePairs, function(pair)
  {
    isSmaller <- trainingMatrix[secondClass, S4Vectors::first(pair)] < trainingMatrix[secondClass, S4Vectors::second(pair)]
    if(sum(isSmaller) < length(isSmaller) / 2)
      pair[2:1]
    else
      pair
  }))
  
  weightingText <- weighted
  if(weightingText == "both") weightingText <- c("unweighted", "weighted")
  testDataFrame <- data.frame(t(testingMatrix), check.names = FALSE)
  
  if(verbose == 3)
    message("Predicting sample classes using feature pair inequalities.")
  
  testPredictions <- do.call(rbind, lapply(testDataFrame, function(sampleValues)
  {
    names(sampleValues) <- rownames(testDataFrame)
    do.call(rbind, lapply(weightingText, function(isWeighted)
    {
      do.call(rbind, lapply(minDifference, function(difference)
      {
        measureDifferences <- sampleValues[S4Vectors::second(featurePairs)] - sampleValues[S4Vectors::first(featurePairs)]
        useFeatures <- which(abs(measureDifferences) > difference)
        if(length(useFeatures) == 0) # No features have a large enough distance difference.
        {                            # Simply vote for the larger class.
          if(largerClass == levels(classes)[1])
          {
            class <- levels(classes)[1]
            score <- -1
          } else {
            class <- levels(classes)[2]
            score <- 1
          }
        } else { # One or more features are available to vote with.
          measureDifferences <- measureDifferences[useFeatures]
          if(isWeighted == "unweighted")
          {
            # For being in second class.
            score <- sum(measureDifferences > 0)
            
            if(score > length(measureDifferences) / 2)
              class <- levels(classes)[2]
            else
              class <- levels(classes)[1]
            
          } else { # Each pair contributes a score for class prediction.
            # For being in second class.
            score <- sum(measureDifferences)

            # Sum of scores is tested for being positive or negative.
            class <- levels(classes)[(sum(measureDifferences) > 0) + 1]
          }
        }
          data.frame(class = factor(class, levels = levels(classes)), score = score,
                     weighted = isWeighted, minDifference = difference)
      }))
    }))
  }))
  
  whichVarieties <- character()
  if(weighted == "both") whichVarieties <- "weighted"
  if(length(minDifference) > 1) whichVarieties <- c(whichVarieties, "minDifference")
  if(length(whichVarieties) == 0) whichVarieties <- "minDifference" # Aribtrary, to make a list.

  varietyFactor <- do.call(paste, c(lapply(whichVarieties, function(variety) paste(variety, testPredictions[, variety], sep = '=')), sep = ','))
  varietyFactor <- factor(varietyFactor)
  resultsList <- lapply(levels(varietyFactor), function(variety)
  {
    varietyPredictions <- subset(testPredictions, varietyFactor == variety)
    
    switch(returnType, class = varietyPredictions[, "class"],
           score = varietyPredictions[, "score"],
           both = data.frame(class = varietyPredictions[, "class"], score = varietyPredictions[, "score"]))
  })
  names(resultsList) <- levels(varietyFactor)
  
  if(length(resultsList) == 1) # No varieties.
    resultsList[[1]]
  else
    resultsList  
})

setMethod("kTSPclassifier", "MultiAssayExperiment", 
          function(measurements, test, target = names(measurements)[1], featurePairs = NULL, ...)
{
  if(is.null(featurePairs))
    stop("No feature pairs provided but some must be.")
  if(!"Pairs" %in% class(featurePairs))
    stop("'featurePairs' must be of type Pairs.")            
            
  tablesAndClasses <- .MAEtoWideTable(measurements, target)
  trainingMatrix <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  testingMatrix <- .MAEtoWideTable(test, target)
            
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  kTSPclassifier(trainingMatrix, classes, testingMatrix, featurePairs, ...)
})