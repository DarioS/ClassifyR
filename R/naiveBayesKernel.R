setGeneric("naiveBayesKernel", function(expression, ...)
           {standardGeneric("naiveBayesKernel")})

setMethod("naiveBayesKernel", "matrix", 
          function(expression, classes, ...)
{ 
  colnames(expression) <- NULL # Might be duplicates because of sampling with replacement.            
  features <- rownames(expression)
  groupsTable <- data.frame(class = classes)  
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features
  naiveBayesKernel(exprSet, ...)
})

setMethod("naiveBayesKernel", "ExpressionSet", 
          function(expression, test, densityFunction = density,
                   densityParameters = list(bw = "nrd0", n = 1024, from = expression(min(featureValues)),
                                                                 to = expression(max(featureValues))),
                   weighted = c("both", "unweighted", "weighted"),
                   weight = c("all", "height difference", "crossover distance", "sum differences"),
                   minDifference = 0, returnType = c("label", "score", "both"), verbose = 3)
{
  weighted <- match.arg(weighted)
  weight <- match.arg(weight)
  returnType <- match.arg(returnType)
            
  classes <- pData(expression)[, "class"]
  classesSizes <- sapply(levels(classes), function(class) sum(classes == class))
  largerClass <- names(classesSizes)[which.max(classesSizes)]
  expression <- exprs(expression)      
  if(class(test) == "ExpressionSet") test <- exprs(test)
  
  if(verbose == 3)
    message("Fitting densities.")
  
  densities <- apply(expression, 1, function(featureValues) 
  {
    oneClassExpression <- featureValues[classes == levels(classes)[1]]
    otherClassExpression <- featureValues[classes == levels(classes)[2]]
    densityParameters <- lapply(densityParameters, function(parameter) eval(parameter))
    oneDensity <- do.call(densityFunction, c(list(oneClassExpression), densityParameters))
    otherDensity <- do.call(densityFunction, c(list(otherClassExpression), densityParameters))
    
    list(oneClass = oneDensity, otherClass = otherDensity)
  })
  
  splines <- lapply(densities, function(featureDensities) list(oneClass = splinefun(featureDensities[["oneClass"]][['x']], featureDensities[["oneClass"]][['y']], "natural"),
                                                               otherClass = splinefun(featureDensities[["otherClass"]][['x']], featureDensities[["otherClass"]][['y']], "natural")))

  posteriorsList <- list()
  if(weight != "height difference" && weighted != "unweighted") # Calculate the crossover distance.
  {
    if(verbose == 3)
      message("Calculating crossover points of class densities.")
    
    # Score for second class level.
    posteriorsHorizontal <- mapply(function(featureDensities, featureSplines, testSamples)
    {
      crosses <- .densityCrossover(featureDensities[["oneClass"]], featureDensities[["otherClass"]])
      otherScores <- classesSizes[2] * featureSplines[["otherClass"]](testSamples) - classesSizes[1] * featureSplines[["oneClass"]](testSamples)
      classPredictions <- ifelse(otherScores > 0, levels(classes)[2], levels(classes)[1])
      classScores <- sapply(testSamples, function(testSample) min(abs(testSample - crosses)))
      classScores <- mapply(function(score, prediction) if(prediction == levels(classes)[1]) -score else score, classScores, classPredictions)
      classScores
    }, densities, splines, as.data.frame(t(test)))

    if(weight != "sum differences")
    {
      posteriorsList[[1]] <- t(posteriorsHorizontal)
      names(posteriorsList) = "crossover distance"
    }
  }
  
  if(weight != "crossover distance") # Calculate the height difference.
  {
    if(verbose == 3)
      message("Calculating vertical differences between densities.")
    
    posteriorsVertical <- mapply(function(featureDensities, featureSplines, testSamples)
    {
      classesSizes[2] * featureSplines[["otherClass"]](testSamples) - classesSizes[1] * featureSplines[["oneClass"]](testSamples)
    }, densities, splines, as.data.frame(t(test)))
    
    if(weight != "sum differences")
      posteriorsList <- c(posteriorsList, `height difference` = list(t(posteriorsVertical)))
  }
  
  if(weight %in% c("sum differences", "all") && weighted != "unweighted") # Sum of the horizontal and vertical distances.
  {
    posteriorsList <- c(posteriorsList, `sum differences` = list(t(posteriorsHorizontal) + t(posteriorsVertical)))
  }

  if(verbose == 3)
  {
    switch(returnType, label = ,
                       both = message("Calculating class scores and determining class labels."),
                       score = message("Calculating class scores.")
          )
  }
  
  weightingText <- weighted
  if(weightingText == "both") weightingText <- c("unweighted", "weighted")
  testPredictions <- do.call(rbind, mapply(function(weightPredictions, weightNames)
  {
    do.call(rbind, lapply(weightingText, function(isWeighted)
    {
      do.call(rbind, lapply(minDifference, function(difference)
      {
        do.call(rbind, apply(weightPredictions, 2, function(sampleCol)
        {
          sampleCol <- sampleCol[abs(sampleCol) > difference]
          if(length(sampleCol) == 0) # No features have a large enough density difference.
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
            if(isWeighted == "unweighted")
            {
              # For being in second class.
              class <- levels(classes)[(sum(sampleCol > 0) > length(sampleCol) / 2) + 1]
              score <- sum(sampleCol > 0) / length(sampleCol)
            } else {
              # For being in second class.
              class <- levels(classes)[(sum(sampleCol) > 0) + 1]
              score <- sum(sampleCol)
            }
          }
          data.frame(class = factor(class, levels = levels(classes)), score = score,
                     weighted = isWeighted, weight = weightNames,
                     minDifference = difference)
        }))
      }))
    }))
  }, posteriorsList, names(posteriorsList), SIMPLIFY = FALSE))
  
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
                 })

  attr(resultsList, "class") <- "list"
  attr(resultsList, "call") <- NULL

  if(length(resultsList) == 1) # No varieties.
    resultsList[[1]]
  else
    resultsList
})