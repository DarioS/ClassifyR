setGeneric("naiveBayesKernel", function(measurements, ...)
           {standardGeneric("naiveBayesKernel")})

setMethod("naiveBayesKernel", "matrix", 
          function(measurements, classes, test, ...)
{
.naiveBayesKernel(DataFrame(t(measurements[, , drop = FALSE]), check.names = FALSE),
                  classes,
                  DataFrame(t(test[, , drop = FALSE]), check.names = FALSE), ...)
})

setMethod("naiveBayesKernel", "DataFrame", 
          function(measurements, classes, test, ...)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  trainingMatrix <- as.matrix(splitDataset[["measurements"]])
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  isNumeric <- sapply(test, is.numeric)
  testingMatrix <- as.matrix(test[, isNumeric, drop = FALSE])
  
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  .naiveBayesKernel(trainingMatrix, splitDataset[["classes"]], testingMatrix, ...)
})

setMethod("naiveBayesKernel", "MultiAssayExperiment", 
          function(measurements, test, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  trainingMatrix <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  testingMatrix <- .MAEtoWideTable(test, targets)
            
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  .naiveBayesKernel(trainingMatrix, classes, testingMatrix, ...)
})

.naiveBayesKernel <- function(measurements, classes, test, densityFunction = density,
                              densityParameters = list(bw = "nrd0", n = 1024, from = expression(min(featureValues)),
                                                       to = expression(max(featureValues))),
                              weighted = c("both", "unweighted", "weighted"),
                              weight = c("all", "height difference", "crossover distance", "sum differences"),
                              minDifference = 0, returnType = c("label", "score", "both"), verbose = 3)
{
  weighted <- match.arg(weighted)
  weight <- match.arg(weight)
  returnType <- match.arg(returnType)
  
  classesSizes <- sapply(levels(classes), function(class) sum(classes == class))
  largerClass <- names(classesSizes)[which.max(classesSizes)]
  
  if(verbose == 3)
    message("Fitting densities.")

  densities <- apply(measurements, 2, function(featureValues)
  {
    oneClassMeasurements <- featureValues[classes == levels(classes)[1]]
    otherClassMeasurements <- featureValues[classes == levels(classes)[2]]
    densityParameters <- lapply(densityParameters, function(parameter) eval(parameter))
    oneDensity <- do.call(densityFunction, c(list(oneClassMeasurements), densityParameters))
    otherDensity <- do.call(densityFunction, c(list(otherClassMeasurements), densityParameters))
    
    list(oneClass = oneDensity, otherClass = otherDensity)
  })
  
  splines <- lapply(densities, function(featureDensities) list(oneClass = splinefun(featureDensities[["oneClass"]][['x']], featureDensities[["oneClass"]][['y']], "natural"),
                                                               otherClass = splinefun(featureDensities[["otherClass"]][['x']], featureDensities[["otherClass"]][['y']], "natural")))
  
  if(weight != "height difference" && weighted != "unweighted") # Calculate the crossover distance.
  {
    if(verbose == 3)
      message("Calculating horizontal distances to crossover points of class densities.")
    
    # Score for second class level.
    distancesHorizontal <- t(do.call(cbind, mapply(function(featureDensities, featureSplines, testSamples)
    {
      crosses <- .densityCrossover(featureDensities[["oneClass"]], featureDensities[["otherClass"]])
      otherScores <- classesSizes[2] * featureSplines[["otherClass"]](testSamples) - classesSizes[1] * featureSplines[["oneClass"]](testSamples)
      classPredictions <- ifelse(otherScores > 0, levels(classes)[2], levels(classes)[1])
      classScores <- sapply(testSamples, function(testSample) min(abs(testSample - crosses)))
      classScores <- mapply(function(score, prediction) if(prediction == levels(classes)[1]) -score else score, classScores, classPredictions)
      classScores
    }, densities, splines, as.data.frame(test), SIMPLIFY = FALSE)))

    class1RelativeSize <- classesSizes[1] / classesSizes[2]
    posteriorsHorizontal <- ifelse(distancesHorizontal < 0, distancesHorizontal * class1RelativeSize, distancesHorizontal * 1/class1RelativeSize)
  }
  
  if(weight != "crossover distance") # Calculate the height difference.
  {
    if(verbose == 3)
      message("Calculating class densities.")
    
    # Rows are for features, columns are for samples.
    distancesVertical <- t(do.call(cbind, mapply(function(featureDensities, featureSplines, testSamples)
    {
      featureSplines[["otherClass"]](testSamples) - featureSplines[["oneClass"]](testSamples)
    }, densities, splines, as.data.frame(test), SIMPLIFY = FALSE)))
    
    # Rows are for features, columns are for samples.
    posteriorsVertical <- t(do.call(cbind, mapply(function(featureDensities, featureSplines, testSamples)
    {
      classesSizes[2] * featureSplines[["otherClass"]](testSamples) - classesSizes[1] * featureSplines[["oneClass"]](testSamples)
    }, densities, splines, as.data.frame(test), SIMPLIFY = FALSE)))
  }

  if(verbose == 3)
  {
    switch(returnType, label = ,
           both = message("Calculating class scores and determining class labels."),
           score = message("Calculating class scores.")
    )
  }
  
  if(weight == "all")
    weightExpanded <-  c("height difference", "crossover distance", "sum differences")
  else weightExpanded <- weight
  allPosteriors <- lapply(weightExpanded, function(type)
  {
    switch(type, `height difference` = posteriorsVertical,
           `crossover distance` = posteriorsHorizontal,
           `sum differences` = posteriorsVertical + posteriorsHorizontal
    )
  })
  allDistances <- lapply(weightExpanded, function(type)
  {
    switch(type, `height difference` = distancesVertical,
           `crossover distance` = distancesHorizontal,
           `sum differences` = distancesVertical + distancesHorizontal
    )
  })
  
  weightingText <- weighted
  if(weightingText == "both") weightingText <- c("unweighted", "weighted")
  testPredictions <- do.call(rbind, mapply(function(weightPredictions, weightNames, distances)
  {
    do.call(rbind, lapply(weightingText, function(isWeighted)
    {
      do.call(rbind, lapply(minDifference, function(difference)
      {
        do.call(rbind, lapply(1:ncol(weightPredictions), function(sampleCol)
        {
          useFeatures <- abs(distances[, sampleCol]) > difference
          if(length(useFeatures) == 0) # No features have a large enough density difference.
          {                          # Simply vote for the larger class.
            if(largerClass == levels(classes)[1])
            {
              class <- levels(classes)[1]
              score <- -1
            } else {
              class <- levels(classes)[2]
              score <- 1
            }
          } else { # One or more features are available to vote with.
            posteriorsUsed <- weightPredictions[useFeatures, sampleCol]
            if(isWeighted == "unweighted")
            {
              # For being in second class.
              class <- levels(classes)[(sum(posteriorsUsed > 0) > length(posteriorsUsed) / 2) + 1]
              score <- sum(posteriorsUsed > 0) / length(posteriorsUsed)
            } else {
              # For being in second class.
              class <- levels(classes)[(sum(posteriorsUsed) > 0) + 1]
              score <- sum(posteriorsUsed)
            }
          }
          data.frame(class = factor(class, levels = levels(classes)), score = score,
                     weighted = isWeighted, weight = weightNames,
                     minDifference = difference)
        }))
      }))
    }))
  }, allPosteriors, weightExpanded, allDistances, SIMPLIFY = FALSE))
  
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
  if(weight == "all") whichVarieties <- c(whichVarieties, "weight")
  if(length(minDifference) > 1) whichVarieties <- c(whichVarieties, "minDifference")
  if(length(whichVarieties) == 0) whichVarieties <- "minDifference" # Aribtrary, to make a list.
  
  varietyFactor <- factor(do.call(paste, c(lapply(whichVarieties, function(variety) paste(variety, testPredictions[, variety], sep = '=')), sep = ',')))
  varietyFactor <- gsub("(weighted=unweighted),weight=height difference", "\\1", varietyFactor)
  resultsList <- by(testPredictions, varietyFactor, function(predictionSet)
  {
    switch(returnType, label = predictionSet[, "class"],
           score = predictionSet[, "score"],
           both = data.frame(label = predictionSet[, "class"], score = predictionSet[, "score"]))
  }, simplify = FALSE)
  
  attr(resultsList, "class") <- "list"
  attr(resultsList, "call") <- NULL
  
  if(length(resultsList) == 1) # No varieties.
    resultsList[[1]]
  else
    resultsList  
}