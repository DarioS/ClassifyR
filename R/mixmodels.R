setGeneric("mixModelsTrain", function(measurements, ...)
           standardGeneric("mixModelsTrain"))

setMethod("mixModelsTrain", "matrix", # Matrix of numeric measurements.
          function(measurements, ...)
{
  mixModelsTrain(DataFrame(t(measurements), check.names = FALSE), ...)
})

setMethod("mixModelsTrain", "DataFrame", # Mixed data types.
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
  
  models <- lapply(levels(classes), function(class)
            {
                aClassMeasurements <- measurements[classes == class, , drop = FALSE]
                apply(aClassMeasurements, 2, function(featureColumn)
                {
                   mixmodParams <- list(featureColumn)
                   mixmodParams <- append(mixmodParams, list(...))
                   do.call(Rmixmod::mixmodCluster, mixmodParams)
                })     
            })
  
  if(verbose == 3)
    message("Done fitting normal mixtures.")
  
  models <- lapply(models, function(modelSet)
            {
              class(modelSet) <- "MixModelsList"
              modelSet
            })
  names(models)[1:length(levels(classes))] <- paste(levels(classes), "Models", sep = '')
  models[["classSizes"]] <- setNames(as.vector(table(classes)), levels(classes))
  models <- MixModelsListsSet(models)
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
           standardGeneric("mixModelsPredict"))

setMethod("mixModelsPredict", c("MixModelsListsSet", "matrix"), function(models, test, ...)
{
  mixModelsPredict(models, DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("mixModelsPredict", c("MixModelsListsSet", "DataFrame"), # Clinical data only.
          function(models, test, weighted = c("unweighted", "weighted", "both"),
                   weight = c("height difference", "crossover distance", "both"),
                   densityXvalues = 1024, minDifference = 0,
                   returnType = c("class", "score", "both"), verbose = 3)
{
  models <- models@set
  isNumeric <- sapply(test, is.numeric)
  test <- test[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  
  weighted <- match.arg(weighted)
  weight <- match.arg(weight)
  returnType <- match.arg(returnType)
  classesNames <- names(models[["classSizes"]])
  classesSizes <- models[["classSizes"]]
  largestClass <- names(classesSizes)[which.max(classesSizes)[1]]
  models <- models[-length(models)]

  if(verbose == 3)
    message("Predicting using normal mixtures.")

  featuresDensities <- lapply(1:ncol(test), function(featureIndex)
  {
    featureValues <- unlist(lapply(models, function(classModels) classModels[[featureIndex]]@data))
    xValues <- seq(min(featureValues), max(featureValues), length.out = densityXvalues)
    setNames(lapply(models, function(model)
    {
      yValues <- Reduce('+', lapply(1:model[[featureIndex]]@bestResult@nbCluster, function(index)
      {
        model[[featureIndex]]@bestResult@parameters@proportions[index] * dnorm(xValues, model[[featureIndex]]@bestResult@parameters@mean[index], sqrt(as.numeric(model[[featureIndex]]@bestResult@parameters@variance[[index]])))
      }))
      list(x = xValues, y = yValues)
    }), classesNames)
  })
      
  splines <- lapply(featuresDensities, function(featureDensities)
             {
               lapply(featureDensities, function(classDensities)
               {
                 splinefun(classDensities[['x']], classDensities[['y']], "natural")   
               })
             })
  
  if(verbose == 3)
    message("Calculating vertical differences between normal mixture densities.")

  # Needed even if horizontal distance weighting is used to determine the predicted class.
  posteriorsVertical <- mapply(function(featureSplines, testSamples)
  {
    sapply(1:length(classesNames), function(classIndex)
    {
      featureSplines[[classIndex]](testSamples)
    })
  }, splines, test, SIMPLIFY = FALSE)
    
  classesVertical <- sapply(posteriorsVertical, function(featureVertical)
  {
      apply(featureVertical, 1, function(sampleVertical) classesNames[which.max(sampleVertical)])
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
    
    classesVerticalIndices <- matrix(match(classesVertical, classesNames),
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
                       score = message("Calculating class scores.")
    )
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
            classScores <- classesSizes / sum(classesSizes)
          } else { # One or more features are available to vote with.
            distancesUsed <- distances[sampleRow, useFeatures]
            classPredictionsUsed <- factor(classesVertical[sampleRow, useFeatures], classesNames)
            if(isWeighted == "unweighted")
            {
              classScores <- table(classPredictionsUsed)
              classScores <- setNames(as.vector(classScores), classesNames)
            } else { # Weighted voting.
              classScores <- tapply(distancesUsed, classPredictionsUsed, sum)
              classScores[is.na(classScores)] <- 0
            }
            classScores <- classScores / sum(classScores) # Make different feature selection sizes comparable.
            classPredicted <- names(classScores)[which.max(classScores)]
          }

          data.frame(class = factor(classPredicted, levels = classesNames), t(classScores),
                     weighted = isWeighted, weight = weightNames,
                     minDifference = difference)
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
  if(length(whichVarieties) == 0) whichVarieties <- "minDifference" # Aribtrary, to make a list.
  
  varietyFactor <- do.call(paste, c(lapply(whichVarieties, function(variety) paste(variety, testPredictions[, variety], sep = '=')), sep = ','))
  varietyFactor <- factor(gsub("(weighted=unweighted),weight=height difference", "\\1", varietyFactor))
  resultsList <- lapply(levels(varietyFactor), function(variety)
  {
    varietyPredictions <- subset(testPredictions, varietyFactor == variety)
    rownames(varietyPredictions) <- rownames(test)
    switch(returnType, class = varietyPredictions[, "class"],
           score = varietyPredictions[, colnames(varietyPredictions) %in% classesNames],
           both = data.frame(class = varietyPredictions[, "class"], varietyPredictions[, colnames(varietyPredictions) %in% classesNames])
           )
  })
  names(resultsList) <- levels(varietyFactor)
  
  if(length(resultsList) == 1) # No varieties.
    resultsList[[1]]
  else
    resultsList  
})

# One or more omics data sets, possibly with clinical data.
setMethod("mixModelsPredict", c("MixModelsListsSet", "MultiAssayExperiment"),
          function(models, test, targets = names(test), ...)
{
  testingMatrix <- .MAEtoWideTable(test, targets)
  mixModelsPredict(models, testingMatrix, ...)
})
