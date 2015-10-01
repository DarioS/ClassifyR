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
                   densityParameters = list(bw = "SJ", n = 1024, from = expression(min(featureValues)),
                                                                 to = expression(max(featureValues))),
                   weighted = c("both", "unweighted", "weighted"),
                   weight = c("both", "height difference", "crossover distance"),
                   minDifference = 0, tolerance = 0.01, returnType = c("label", "score", "both"), verbose = 3)
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
  if(weight != "height difference") # Calculate the crossover distance.
  {
    if(verbose == 3)
      message("Calculating crossover points of class densities.")
    
    # Score for second class level.
    posteriorsHorizontal <- mapply(function(featureDensities, featureSplines, testSamples)
    {
      crosses <- .densityCrossover(featureDensities[["oneClass"]], featureDensities[["otherClass"]], tolerance)
      otherScores <- classesSizes[2] * featureSplines[["otherClass"]](testSamples) - classesSizes[1] * featureSplines[["oneClass"]](testSamples)
      classPredictions <- ifelse(otherScores > 0, levels(classes)[2], levels(classes)[1])
      classScores <- sapply(testSamples, function(testSample) min(abs(testSample - crosses)))
      classScores <- mapply(function(score, prediction) if(prediction == levels(classes)[1]) -score else score, classScores, classPredictions)
      classScores
    }, densities[1], splines[1], as.data.frame(t(test))[1])

    posteriorsList[[1]] <- t(posteriorsHorizontal)
    names(posteriorsList) = "crossover distance"
  }
  
  if(weight != "crossover distance") # Calculate the height difference.
  {
    if(verbose == 3)
      message("Calculating vertical differences between densities.")
    
    posteriorsVertical <- mapply(function(featureDensities, featureSplines, testSamples)
    {
      classesSizes[2] * featureSplines[["otherClass"]](testSamples) - classesSizes[1] * featureSplines[["oneClass"]](testSamples)
    }, densities, splines, as.data.frame(t(test)))
    
    posteriorsList <- c(posteriorsList, `height difference` = list(t(posteriorsVertical)))
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
              logicalSymbol <- FALSE
              score <- -1
            } else {
              logicalSymbol <- TRUE
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
          data.frame(class = class, score = score,
                     weighted = isWeighted, weight = weightNames,
                     minDifference = difference)
        }))
      }))
    }))
  }, posteriorsList, names(posteriorsList), SIMPLIFY = FALSE))
  
  whichVarieties <- character()
  if(weighted == "both") whichVarieties <- "weighted"
  if(weight == "both") whichVarieties <- c(whichVarieties, "weight")
  if(length(minDifference) > 1) whichVarieties <- c(whichVarieties, "minDifference")
  if(length(whichVarieties) == 0) whichVarieties <- "minDifference" # Aribtrary, to make a list.

  varietyFactor <- factor(do.call(paste, c(lapply(whichVarieties, function(variety) paste(variety, testPredictions[, variety], sep = '=')), sep = ',')))
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