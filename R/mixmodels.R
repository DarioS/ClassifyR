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
          function(models, test, difference = c("unweighted", "weighted"),
                   weighting = c("height difference", "crossover distance"),
                   densityXvalues = 1024, minDifference = 0,
                   returnType = c("both", "class", "score"), verbose = 3)
{
  models <- models@set
  isNumeric <- sapply(test, is.numeric)
  test <- test[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  
  difference <- match.arg(difference)
  weighting <- match.arg(weighting)
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
  
  if(difference == "crossover distance")
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
  
  allDistances <- switch(weighting, `height difference` = distancesVertical,
                         `crossover distance` = distancesHorizontal)
  
  predictions <- do.call(rbind, lapply(1:nrow(allDistances), function(sampleRow)
  {
    useFeatures <- abs(allDistances[sampleRow, ]) > minDifference
    if(all(useFeatures == FALSE)) # No features have a large enough density difference.
    {                          # Simply vote for the larger class.
      classPredicted <- largestClass
      classScores <- classesSizes / length(classes)
    } else { # One or more features are available to vote with.
      distancesUsed <- allDistances[sampleRow, useFeatures]
      classPredictionsUsed <- factor(classesVertical[sampleRow, useFeatures], levels(classes))
      if(difference == "unweighted")
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
    
    data.frame(class = factor(classPredicted, levels = levels(classes)), t(classScores), check.names = FALSE)
  }))
  
  switch(returnType, class = predictions[, "class"],
         score = predictions[, colnames(predictions) %in% levels(classes)],
         both = data.frame(class = predictions[, "class"], predictions[, colnames(predictions) %in% levels(classes)], check.names = FALSE)
  ) 
})

# One or more omics data sets, possibly with clinical data.
setMethod("mixModelsPredict", c("MixModelsListsSet", "MultiAssayExperiment"),
          function(models, test, targets = names(test), ...)
{
  testingMatrix <- .MAEtoWideTable(test, targets)
  mixModelsPredict(models, testingMatrix, ...)
})
