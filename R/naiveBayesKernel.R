setGeneric("naiveBayesKernel", function(measurements, ...)
           standardGeneric("naiveBayesKernel"))

setMethod("naiveBayesKernel", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, test, ...)
{
  naiveBayesKernel(DataFrame(t(measurements[, , drop = FALSE]), check.names = FALSE),
                   classes,
                   DataFrame(t(test[, , drop = FALSE]), check.names = FALSE), ...)
})

setMethod("naiveBayesKernel", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, test,
                   densityFunction = density, densityParameters = list(bw = "nrd0", n = 1024, from = expression(min(featureValues)), to = expression(max(featureValues))),
                   weighted = c("unweighted", "weighted", "both"),
                   weight = c("height difference", "crossover distance", "both"),
                   minDifference = 0, returnType = c("both", "class", "score"), verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  trainingMatrix <- splitDataset[["measurements"]]
  isNumeric <- sapply(trainingMatrix, is.numeric)
  trainingMatrix <- as.matrix(trainingMatrix[, isNumeric, drop = FALSE])
  classes <- splitDataset[["classes"]]
  isNumeric <- sapply(test, is.numeric)
  testingMatrix <- as.matrix(test[, isNumeric, drop = FALSE])
  
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  
  weighted <- match.arg(weighted)
  weight <- match.arg(weight)
  returnType <- match.arg(returnType)
  
  classesSizes <- sapply(levels(classes), function(class) sum(classes == class))
  largestClass <- names(classesSizes)[which.max(classesSizes)[1]]
  
  if(verbose == 3)
    message("Fitting densities.")

  featuresDensities <- lapply(measurements, function(featureValues)
  {
    densityParameters <- lapply(densityParameters, function(parameter) eval(parameter))
    lapply(levels(classes), function(class)
    {
      aClassMeasurements <- featureValues[classes == class]  
      do.call(densityFunction, c(list(aClassMeasurements), densityParameters))
    }) # A fitted density for each class.
  })

  classesScaleFactors <- classesSizes / nrow(measurements)
  splines <- lapply(featuresDensities, function(featureDensities) 
             {
               mapply(function(featureDensity, scaleFactor)
               {
                 splinefun(featureDensity[['x']], featureDensity[['y']] * scaleFactor, "natural")
               }, featureDensities, classesScaleFactors)
             })
  
  if(verbose == 3)
    message("Calculating vertical distances between class densities.")

  # Needed even if horizontal distance weighting is used to determine the predicted class.
  posteriorsVertical <- mapply(function(featureSplines, testSamples)
  {
    sapply(1:length(levels(classes)), function(classIndex)
    {
      featureSplines[[classIndex]](testSamples)
    })
  }, splines, test, SIMPLIFY = FALSE)
    
  classesVertical <- sapply(posteriorsVertical, function(featureVertical)
  {
      apply(featureVertical, 1, function(sampleVertical) levels(classes)[which.max(sampleVertical)])
  }) # Matrix, rows are test samples, columns are features.
    
  distancesVertical <- sapply(posteriorsVertical, function(featureVertical)
  { # Vertical distance between highest density and second-highest, at a particular value.
    apply(featureVertical, 1, function(sampleVertical)
    {
      twoHighest <- sort(sampleVertical, decreasing = TRUE)[1:2]
      Reduce('-', twoHighest)
    })
  }) # Matrix, rows are test samples, columns are features.
  
  if(weight %in% c("crossover distance", "both")) # Calculate the crossover distance, even if unweighted voting to pick the class.
  {
    if(verbose == 3)
      message("Calculating horizontal distances to crossover points of class densities.")
 
    classesVerticalIndices <- matrix(match(classesVertical, levels(classes)),
                                     nrow = nrow(classesVertical), ncol = ncol(classesVertical))
    distancesHorizontal <- mapply(function(featureDensities, testSamples, predictedClasses)
    {
      classesCrosses <- .densitiesCrossover(featureDensities)
      classesDistances <- sapply(classesCrosses, function(classCrosses)
      {
        sapply(testSamples, function(testSample) min(abs(testSample - classCrosses)))
      })
      classesDistances[cbind(1:nrow(classesDistances), predictedClasses)]
    }, featuresDensities, test, as.data.frame(classesVerticalIndices)) # Matrix of horizontal distances to nearest cross-over involving the predicted class.
  }

  if(verbose == 3)
  {
    switch(returnType, class = message("Determining class labels."),
                       both = message("Calculating class scores and determining class labels."),
                       score = message("Calculating class scores."))
  }

  if(weight == "both")
    weightExpanded <-  c("height difference", "crossover distance")
  else weightExpanded <- weight
  
  allDistances <- lapply(weightExpanded, function(type)
  {
    switch(type, `height difference` = distancesVertical,
                 `crossover distance` = distancesHorizontal)
  })

  weightingText <- weighted
  if(weightingText == "both") weightingText <- c("unweighted", "weighted")
  testPredictions <- do.call(rbind, mapply(function(weightNames, distances)
  {
    do.call(rbind, lapply(weightingText, function(isWeighted)
    {
      do.call(rbind, lapply(minDifference, function(difference)
      {
        do.call(rbind, lapply(1:nrow(distances), function(sampleRow)
        {
          useFeatures <- abs(distances[sampleRow, ]) > difference
          if(all(useFeatures == FALSE)) # No features have a large enough density difference.
          {                          # Simply vote for the larger class.
            classPredicted <- largestClass
            classScores <- classesSizes / length(classes)
          } else { # One or more features are available to vote with.
            distancesUsed <- distances[sampleRow, useFeatures]
            classPredictionsUsed <- factor(classesVertical[sampleRow, useFeatures], levels(classes))
            if(isWeighted == "unweighted")
            {
              classScores <- table(classPredictionsUsed)
              classScores <- setNames(as.vector(classScores), levels(classes))
            } else { # Weighted voting.
              classScores <- tapply(distancesUsed, classPredictionsUsed, sum)
              classScores[is.na(classScores)] <- 0
            }
            classScores <- classScores / sum(classScores) # Make different feature selection sizes comparable.
            classPredicted <- names(classScores)[which.max(classScores)]
          }

          data.frame(class = factor(classPredicted, levels = levels(classes)), t(classScores),
                     weighted = isWeighted, weight = weightNames,
                     minDifference = difference, check.names = FALSE)
        }))
      }))
    }))
  }, weightExpanded, allDistances, SIMPLIFY = FALSE))

  # Remove combinations of unweighted voting and weightings.
  testPredictions <- do.call(rbind, by(testPredictions, testPredictions[, "weighted"], function(weightVariety)
  {
    if(weightVariety[1, "weighted"] == "unweighted")
    {
      do.call(rbind, by(weightVariety, weightVariety[, "minDifference"], function(differenceVariety) differenceVariety[differenceVariety[, "weight"] == "height difference", ]))
    } else {
      weightVariety
    }
  }))

  whichVarieties <- character()
  if(weighted == "both") whichVarieties <- "weighted"
  if(weight == "both") whichVarieties <- c(whichVarieties, "weight")
  if(length(minDifference) > 1) whichVarieties <- c(whichVarieties, "minDifference")
  if(length(whichVarieties) == 0) whichVarieties <- "minDifference" # Arbitrary, to make a list.

  varietyFactor <- do.call(paste, c(lapply(whichVarieties, function(variety) paste(variety, testPredictions[, variety], sep = '=')), sep = ','))
  varietyFactor <- factor(gsub("(weighted=unweighted),weight=height difference", "\\1", varietyFactor))
  resultsList <- lapply(levels(varietyFactor), function(variety)
  {
    varietyPredictions <- subset(testPredictions, varietyFactor == variety)
    rownames(varietyPredictions) <- rownames(test)
    switch(returnType, class = varietyPredictions[, "class"],
           score = varietyPredictions[, colnames(varietyPredictions) %in% levels(classes)],
           both = data.frame(class = varietyPredictions[, "class"], varietyPredictions[, colnames(varietyPredictions) %in% levels(classes)], check.names = FALSE)
           )
  })
  names(resultsList) <- levels(varietyFactor)

  if(length(resultsList) == 1) # No varieties.
    resultsList[[1]]
  else
    resultsList    
})

setMethod("naiveBayesKernel", "MultiAssayExperiment", 
          function(measurements, test, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  trainingMatrix <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  testingMatrix <- .MAEtoWideTable(test, targets)
            
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  naiveBayesKernel(trainingMatrix, classes, testingMatrix, ...)
})