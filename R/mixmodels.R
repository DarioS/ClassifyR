setGeneric("mixModelsTrain", function(measurements, ...)
           {standardGeneric("mixModelsTrain")})

setMethod("mixModelsTrain", "matrix", # Matrix of numeric measurements.
          function(measurements, ...)
{
  mixModelsTrain(DataFrame(t(measurements), check.names = FALSE), ...)
})

setMethod("mixModelsTrain", "DataFrame", # Clinical data only.
          function(measurements, classes, ..., verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  
  if(verbose == 3)
    message("Fitting mixtures of normals for genes.")
  if(!requireNamespace("Rmixmod", quietly = TRUE))
    stop("The package 'Rmixmod' could not be found. Please install it.")

  oneClass <- classes == levels(classes)[1]
  otherClass <- classes == levels(classes)[2]
  oneClassMeasurements <- measurements[oneClass, ]
  otherClassMeasurements <- measurements[otherClass, ]
  
  models <- lapply(list(oneClassMeasurements, otherClassMeasurements), function(classMeasurements)
            {
                apply(classMeasurements, 2, function(featureColumn)
                {
                   mixmodParams <- list(featureColumn)
                   mixmodParams <- append(mixmodParams, list(...))
                   do.call(Rmixmod::mixmodCluster, mixmodParams)
                })     
            })
  
  if(verbose == 3)
    message("Done fitting normal mixtures.")
  
  models[["classSizes"]] <- setNames(c(sum(oneClass), sum(otherClass)), levels(classes))
  models
})

# One or more omics data sets, possibly with clinical data.
setMethod("mixModelsTrain", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]            
  mixModelsTrain(dataTable, classes, ...)
})

setGeneric("mixModelsPredict", function(models, test, ...)
           {standardGeneric("mixModelsPredict")})

setMethod("mixModelsPredict", c("list", "matrix"), function(models, test, ...)
{
  mixModelsPredict(models, DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("mixModelsPredict", c("list", "DataFrame"), # Clinical data only.
          function(models, test, weighted = c("both", "unweighted", "weighted"),
                   weight = c("all", "height difference", "crossover distance", "sum differences"),
                   densityXvalues = 1024, minDifference = 0,
                   returnType = c("class", "score", "both"), verbose = 3)
{
  isNumeric <- sapply(test, is.numeric)
  test <- test[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")

  weighted <- match.arg(weighted)
  weight <- match.arg(weight)
  returnType <- match.arg(returnType)
  classLevels <- names(models[["classSizes"]])
  largerClass <- classLevels[which.max(models[["classSizes"]])[1]]
  
  if(verbose == 3)
    message("Predicting using normal mixtures.")

  if(weight != "height difference") # Calculate the crossover distance.
  {
    if(verbose == 3)
      message("Calculating crossover points of normal mixture densities.")
    # Convert to a density, so the crossover points can be calculated.
    densities <- mapply(function(oneClassModel, otherClassModel)
    {
      featureValues <- c(oneClassModel@data, otherClassModel@data)
      xValues <- seq(min(featureValues), max(featureValues), length.out = densityXvalues)
      setNames(lapply(list(oneClassModel, otherClassModel), function(model)
      {
        yValues <- Reduce('+', lapply(1:model@bestResult@nbCluster, function(index)
        {
          model@bestResult@parameters@proportions[index] * dnorm(xValues, model@bestResult@parameters@mean[index], sqrt(as.numeric(model@bestResult@parameters@variance[[index]])))
        }))
        list(x = xValues, y = yValues)
      }), c("oneClass", "otherClass"))
    }, models[[1]], models[[2]], SIMPLIFY = FALSE)
    crosses <- lapply(densities, function(densityPair) .densityCrossover(densityPair[[1]], densityPair[[2]]))
    splines <- lapply(densities, function(featureDensities) list(oneClass = splinefun(featureDensities[["oneClass"]][['x']], featureDensities[["oneClass"]][['y']], "natural"),
                                                                 otherClass = splinefun(featureDensities[["otherClass"]][['x']], featureDensities[["otherClass"]][['y']], "natural")))

    distancesHorizontal <- t(do.call(cbind, mapply(function(featureCrosses, testSamples)
    {
      sapply(testSamples, function(testSample) min(abs(testSample - featureCrosses)))
    }, crosses, test, SIMPLIFY = FALSE)))

    posteriorsHorizontal <- t(do.call(cbind, mapply(function(featureCrosses, featureSplines, testSamples)
    {
      classScores <- models[["classSizes"]][2] * featureSplines[["otherClass"]](testSamples) - models[["classSizes"]][1] * featureSplines[["oneClass"]](testSamples)
      classPredictions <- ifelse(classScores > 0, classLevels[2], classLevels[1])
      classScores <- sapply(testSamples, function(testSample) min(abs(testSample - featureCrosses)))
      classScores <- mapply(function(score, prediction) if(prediction == classLevels[1]) -score else score, classScores, classPredictions)
      classScores
    }, crosses, splines, test, SIMPLIFY = FALSE)))
    class1RelativeSize <- models[["classSizes"]][1] / models[["classSizes"]][2]
    posteriorsHorizontal <- if(class1RelativeSize < 0) posteriorsHorizontal <- posteriorsHorizontal * class1RelativeSize else posteriorsHorizontal <- posteriorsHorizontal * 1/class1RelativeSize
  }
  
  if(weight != "crossover distance") # Calculate the height difference.
  {
    if(verbose == 3)
      message("Calculating vertical differences between normal mixture densities.")
    
    # Rows are for features, columns are for samples.
    distancesVertical <- t(do.call(cbind, mapply(function(featureSplines, testSamples)
    {
      featureSplines[["otherClass"]](testSamples) - featureSplines[["oneClass"]](testSamples)
    }, splines, as.data.frame(test), SIMPLIFY = FALSE)))
    
    # Rows are for features, columns are for samples.
    posteriorsVertical <- t(do.call(cbind, mapply(function(featureSplines, testSamples)
    {
      models[["classSizes"]][2] * featureSplines[["otherClass"]](testSamples) - models[["classSizes"]][1] * featureSplines[["oneClass"]](testSamples)
    }, splines, as.data.frame(test), SIMPLIFY = FALSE)))
  }
  
  if(verbose == 3)
  {
    switch(returnType, class = ,
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
          if(sum(useFeatures) == 0) # No features have a large enough density difference.
          {                          # Simply vote for the larger class.
            if(largerClass == classLevels[1])
            {
              class <- classLevels[1]
              score <- -1
            } else {
              class <- classLevels[2]
              score <- 1
            }
          } else { # One or more features are available to vote with.
            sampleValues <- weightPredictions[useFeatures, sampleCol]
            if(isWeighted == "unweighted")
            {
              # For being in second class.
              class <- classLevels[(sum(sampleValues > 0) > length(sampleValues) / 2) + 1]
              score <- sum(sampleValues > 0) / length(sampleValues)
            } else {
              # For being in second class.
              class <- classLevels[(sum(sampleValues) > 0) + 1]
              score <- sum(sampleValues)
            }
          }
          data.frame(class = factor(class, levels = classLevels), score = score,
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
  
  varietyFactor <- do.call(paste, c(lapply(whichVarieties, function(variety) paste(variety, testPredictions[, variety], sep = '=')), sep = ','))
  varietyFactor <- factor(gsub("(weighted=unweighted),weight=height difference", "\\1", varietyFactor))
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

# One or more omics data sets, possibly with clinical data.
setMethod("mixModelsPredict", c("list", "MultiAssayExperiment"),
          function(models, test, targets = names(test), ...)
{
  testingMatrix <- .MAEtoWideTable(test, targets)
  mixModelsPredict(models, testingMatrix, ...)
})