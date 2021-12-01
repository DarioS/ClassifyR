setGeneric("kTSPclassifier", function(measurements, ...)
           standardGeneric("kTSPclassifier"))

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
                   difference = c("unweighted", "weighted"), minDifference = 0,
                   returnType = c("both", "class", "score"), verbose = 3)
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
  
  difference <- match.arg(difference)
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
      Pairs(S4Vectors::second(pair), S4Vectors::first(pair))
    else
      pair
  }))

  testDataFrame <- data.frame(t(testingMatrix), check.names = FALSE)
  
  if(verbose == 3)
    message("Predicting sample classes using feature pair inequalities.")
  
  predictions <- do.call(rbind, lapply(testDataFrame, function(sampleValues)
  {
    names(sampleValues) <- rownames(testDataFrame)
    measureDifferences <- sampleValues[S4Vectors::second(featurePairs)] - sampleValues[S4Vectors::first(featurePairs)]
    useFeatures <- which(abs(measureDifferences) > minDifference)
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
      if(difference == "unweighted")
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
    data.frame(class = factor(class, levels = levels(classes)), score = score, check.names = FALSE)
  }))

  switch(returnType, class = predictions[, "class"],
         score = predictions[, colnames(predictions) %in% levels(classes)],
         both = data.frame(class = predictions[, "class"], predictions[, colnames(predictions) %in% levels(classes)], check.names = FALSE)
  )
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