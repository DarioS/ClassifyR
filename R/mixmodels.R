setGeneric("mixModelsTrain", function(measurements, ...)
           {standardGeneric("mixModelsTrain")})

setMethod("mixModelsTrain", "matrix", # Matrix of numeric measurements.
          function(measurements, ...)
{
  .mixModelsTrain(DataFrame(t(measurements), check.names = FALSE), ...)
})

setMethod("mixModelsTrain", "DataFrame", # Clinical data only.
          function(measurements, classes, ...)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  .mixModelsTrain(measurements, splitDataset[["classes"]], ...)
})

# One or more omics datasets, possibly with clinical data.
setMethod("mixModelsTrain", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]            
  .mixModelsTrain(dataTable, classes, ...)
})

.mixModelsTrain <- function(measurements, classes, ..., verbose = 3)
{
  if(verbose == 3)
    message("Fitting mixtures of normals for genes.")
  if(!requireNamespace("Rmixmod", quietly = TRUE))
    stop("The package 'Rmixmod' could not be found. Please install it.")

  oneClass <- classes == levels(classes)[1]
  otherClass <- classes == levels(classes)[2]
  oneClassMeasurements <- measurements[oneClass, ]
  otherClassMeasurements <- measurements[otherClass, ]
  isExtra <- !missing(...)
  
  models <- lapply(list(oneClassMeasurements, otherClassMeasurements), function(classMeasurements)
            {
                apply(classMeasurements, 2, function(featureColumn)
                {
                   mixmodParams <- list(featureColumn)
                   if(isExtra) mixmodParams <- append(mixmodParams, list(...))
                   do.call(Rmixmod::mixmodCluster, mixmodParams)
                })     
            })
  
  if(verbose == 3)
    message("Done fitting normal mixtures.")
  
  models[["classSizes"]] <- setNames(c(sum(oneClass), sum(otherClass)), levels(classes))
  models
}

setGeneric("mixModelsTest", function(models, test, ...)
           {standardGeneric("mixModelsTest")})

setMethod("mixModelsTest", c("list", "matrix"), function(models, test, ...)
{
  .mixModelsTest(models, DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("mixModelsTest", c("list", "DataFrame"), # Clinical data only.
          function(models, test, ...)
{
  splitDataset <- .splitDataAndClasses(test, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  .mixModelsTest(models, measurements, ...)
})

# One or more omics datasets, possibly with clinical data.
setMethod("mixModelsTest", c("list", "MultiAssayExperiment"),
          function(models, test, targets = names(test), ...)
{
  testingMatrix <- .MAEtoWideTable(test, targets)
  .mixModelsTest(models, testingMatrix, ...)
})

.mixModelsTest <- function(models, test, weighted = c("both", "unweighted", "weighted"),
                  weight = c("all", "height difference", "crossover distance", "sum differences"),
                  densityXvalues = 1024, minDifference = 0,
                  returnType = c("label", "score", "both"), verbose = 3)
{
  weighted <- match.arg(weighted)
  weight <- match.arg(weight)
  returnType <- match.arg(returnType)

  classLevels <- names(models[["classSizes"]])
  
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
    }, models[[1]], models[[2]], SIMPLIFY = FALSE)
    
    crosses <- lapply(densities, function(densityPair) .densityCrossover(densityPair[[1]], densityPair[[2]]))
    
    posteriorsHorizontal <- do.call(cbind, mapply(function(testSamples, oneClassModel, otherClassModel, featureCrosses)
    {
      featureValues <- c(oneClassModel@data, otherClassModel@data)
      classScores <- lapply(list(oneClassModel, otherClassModel), function(model)
      {
        Reduce('+', lapply(1:model@bestResult@nbCluster, function(index)
        {
          model@bestResult@parameters@proportions[index] * dnorm(testSamples, model@bestResult@parameters@mean[index], sqrt(as.numeric(model@bestResult@parameters@variance[[index]])))
        }))
      })
      
      classScores <- models[["classSizes"]][2] * classScores[[2]] - models[["classSizes"]][1] * classScores[[1]]
      classPredictions <- ifelse(classScores > 0, classLevels[2], classLevels[1])
      classScores <- sapply(testSamples, function(testSample) min(abs(testSample - featureCrosses)))
      classScores <- mapply(function(score, prediction) if(prediction == levels(classes)[1]) -score else score, classScores, classPredictions)
      classScores
    }, test, models[[1]], models[[2]], crosses, SIMPLIFY = FALSE))
    
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
    posteriorsVertical <- do.call(cbind, mapply(function(featureValues, oneClassModel, otherClassModel)
                          {
                            classScores <- lapply(list(oneClassModel, otherClassModel), function(model)
                                           {
                                             Reduce('+', lapply(1:model@bestResult@nbCluster, function(index)
                                             {
                                               model@bestResult@parameters@proportions[index] * dnorm(featureValues, model@bestResult@parameters@mean[index], sqrt(as.numeric(model@bestResult@parameters@variance[[index]])))
                                             }))
                                            })
                            classScores <- models[["classSizes"]][2] * classScores[[2]] - models[["classSizes"]][1] * classScores[[1]]
                           }, test, models[[1]], models[[2]], SIMPLIFY = FALSE))
    
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
  }, simplify = FALSE)
  attr(resultsList, "class") <- "list"
  attr(resultsList, "call") <- NULL
  
  if(length(resultsList) == 1) # No varieties.
    resultsList[[1]]
  else
    resultsList
}