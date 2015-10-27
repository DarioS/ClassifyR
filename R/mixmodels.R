setGeneric("mixModelsTrain", function(expression, ...)
           {standardGeneric("mixModelsTrain")})

setMethod("mixModelsTrain", "matrix", function(expression, classes, ...)
{ 
  colnames(expression) <- NULL # Might be duplicates because of sampling with replacement.
  features <- rownames(expression)
  groupsTable <- data.frame(class = classes)
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features
  mixModelsTrain(exprSet, ...)
})

setMethod("mixModelsTrain", "ExpressionSet",
          function(expression, ..., verbose = 3)
{
  if(verbose == 3)
    message("Fitting mixtures of normals for genes.")
  if(!requireNamespace("Rmixmod", quietly = TRUE))
    stop("The package 'Rmixmod' could not be found. Please install it.")
  classes <- pData(expression)[, "class"]
  expression <- exprs(expression)
  oneClass <- classes == levels(classes)[1]
  otherClass <- classes == levels(classes)[2]
  oneClassExpression <- expression[, oneClass]
  otherClassExpression <- expression[, otherClass]
  isExtra <- !missing(...)
  
  models <- lapply(list(oneClassExpression, otherClassExpression), function(classExpression)
            {
                classExpression <- data.frame(t(classExpression))
                list(apply(classExpression, 2, function(geneRow)
                {
                   mixmodParams <- list(geneRow)
                   if(isExtra) mixmodParams <- append(mixmodParams, list(...))
                   do.call(Rmixmod::mixmodCluster, mixmodParams)
                }), nrow(classExpression))        
            })
  
  if(verbose == 3)
    message("Done fitting normal mixtures.")
  attr(models, "classes") <- levels(classes)
  models
})

setGeneric("mixModelsTest", function(models, test, ...)
           {standardGeneric("mixModelsTest")})

setMethod("mixModelsTest", c("list", "matrix"), function(models, test, ...)
{ 
  colnames(test) <- NULL # Might be duplicates because of sampling with replacement.  
  exprSet <- ExpressionSet(test)
  mixModelsTest(models, exprSet, ...)
})

setMethod("mixModelsTest", c("list", "ExpressionSet"),
          function(models, test, weighted = c("both", "unweighted", "weighted"),
                   weight = c("all", "height difference", "crossover distance", "sum differences"),
                   densityXvalues = 1024, minDifference = 0,
                   returnType = c("label", "score", "both"), verbose = 3)
{
  weighted <- match.arg(weighted)
  weight <- match.arg(weight)
  returnType <- match.arg(returnType)
  test <- exprs(test)
  classLevels <- attr(models, "classes")
  
  if(verbose == 3)
    message("Predicting using normal mixtures.")
  
  posteriorsList <- list()
  if(weight != "height difference") # Calculate the crossover distance.
  {
    if(verbose == 3)
      message("Calculating crossover points of normal mixture densities.")
    # Convert to a density, so the crossover points can be calculated.
    densities <- mapply(function(oneClassModel, otherClassModel)
    {
      featureValues <- c(oneClassModel@data, otherClassModel@data)
      xValues <- seq(min(featureValues), max(featureValues), length.out = densityXvalues)
      classScores <- lapply(list(oneClassModel, otherClassModel), function(model)
      {
        yValues <- Reduce('+', lapply(1:model@bestResult@nbCluster, function(index)
        {
          model@bestResult@parameters@proportions[index] * dnorm(xValues, model@bestResult@parameters@mean[index], sqrt(as.numeric(model@bestResult@parameters@variance[[index]])))
        }))
        list(x = xValues, y = yValues)
      })
    }, models[[1]][[1]], models[[2]][[1]], SIMPLIFY = FALSE)
    
    crosses <- lapply(densities, function(densityPair) .densityCrossover(densityPair[[1]], densityPair[[2]]))
    
    posteriorsHorizontal <- mapply(function(testSamples, oneClassModel, otherClassModel, featureCrosses)
    {
      featureValues <- c(oneClassModel@data, otherClassModel@data)
      classScores <- lapply(list(oneClassModel, otherClassModel), function(model)
      {
        Reduce('+', lapply(1:model@bestResult@nbCluster, function(index)
        {
          model@bestResult@parameters@proportions[index] * dnorm(testSamples, model@bestResult@parameters@mean[index], sqrt(as.numeric(model@bestResult@parameters@variance[[index]])))
        }))
      })
      
      classScores <- models[[2]][[2]] * classScores[[2]] - models[[1]][[2]] * classScores[[1]]
      classPredictions <- ifelse(classScores > 0, classLevels[2], classLevels[1])
      classScores <- sapply(testSamples, function(testSample) min(abs(testSample - featureCrosses)))
      classScores <- mapply(function(score, prediction) if(prediction == levels(classes)[1]) -score else score, classScores, classPredictions)
      classScores
      # Second element of second list in 'models' is unimportant information added by mixmod.
    }, data.frame(t(test)), models[[1]][[1]], models[[2]][[1]], crosses)
    
    if(weight != "sum differences")
    {
      posteriorsList[[1]] <- t(posteriorsHorizontal)
      names(posteriorsList) = "crossover distance"
    }
  }
  
  if(weight != "crossover distance") # Calculate the height difference.
  {
    if(verbose == 3)
      message("Calculating vertical differences between normal mixture densities.")
    posteriorsVertical <- mapply(function(featureValues, oneClassModel, otherClassModel)
                          {
                            classScores <- lapply(list(oneClassModel, otherClassModel), function(model)
                                           {
                                             Reduce('+', lapply(1:model@bestResult@nbCluster, function(index)
                                             {
                                               model@bestResult@parameters@proportions[index] * dnorm(featureValues, model@bestResult@parameters@mean[index], sqrt(as.numeric(model@bestResult@parameters@variance[[index]])))
                                             }))
                                            })
                            classScores <- models[[2]][[2]] * classScores[[2]] - models[[1]][[2]] * classScores[[1]]

                            # Second element of second list in 'models' is unimportant information added by mixmod.
                           }, data.frame(t(test)), models[[1]][[1]], models[[2]][[1]])
    
    if(weight != "sum differences")
      posteriorsList <- c(posteriorsList, `height difference` = list(t(posteriorsVertical)))
  }
  
  if(weight == "sum differences") # Sum of the horizontal and vertical distances.
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
