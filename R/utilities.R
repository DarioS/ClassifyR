.splitDataAndClasses <- function(measurements, classes)
{ # DataFrame methods' class variable can be character or factor, so it's a bit involved.
  if(class(classes) == "character" && length(classes) > 1)
    stop("'classes' is a character variable but has more than one element. Either provide a\n",
         "       single column name or a factor of the same length as the number of samples.")
  
  if(class(classes) == "character")
  {
    classColumn <- match(classes, colnames(measurements))
    if(is.na(classColumn))
      stop("Specified column name of classes is not present in the data table.")
    classes <- measurements[, classColumn]
    measurements <- measurements[, -classColumn]
    if(class(classes) != "factor")
      classes <- factor(classes)
  }
  list(measurements = measurements, classes = classes)
}

.MAEtoWideTable <- function(measurements, targets, restrict = "numeric")
{
  if(is.null(targets))
    stop("'targets' is not specified but must be.")  
  if(!all(targets %in% c(names(measurements), "clinical")))
    stop("Some table names in 'targets' are not assay names in 'measurements' or \"clinical\".")  

  if("clinical" %in% targets)
  {
    clinicalColumns <- colnames(MultiAssayExperiment::colData(measurements))
    targets <- targets[-match("clinical", targets)]
  } else if("class" %in% colnames(MultiAssayExperiment::colData(measurements))) {
    clinicalColumns <- "class"
  } else {
    clinicalColumns <- NULL
  }

  if(length(targets) > 0)
  {
    measurements <- measurements[, , targets]
  
    dataTable <- wideFormat(measurements, colDataCols = clinicalColumns, check.names = FALSE, collapse = ':')
    rownames(dataTable) <- dataTable[, "primary"]
    S4Vectors::mcols(dataTable)[, "sourceName"] <- gsub("colDataCols", "clinical", S4Vectors::mcols(dataTable)[, "sourceName"])
    colnames(S4Vectors::mcols(dataTable))[1] <- "dataset"
  
    S4Vectors::mcols(dataTable)[, "feature"] <- as.character(S4Vectors::mcols(dataTable)[, "rowname"])
    missingIndices <- is.na(S4Vectors::mcols(dataTable)[, "feature"])
    S4Vectors::mcols(dataTable)[missingIndices, "feature"] <- colnames(dataTable)[missingIndices]
    S4Vectors::mcols(dataTable) <- S4Vectors::mcols(dataTable)[, c("dataset", "feature")]
    if("class" %in% colnames(dataTable))
      classes <- dataTable[, "class"]
    else
      classes <- NULL
  } else { # Must have only been clinical data.
    dataTable <- MultiAssayExperiment::colData(measurements)
    classes <- dataTable[, "class"]
  }
  
  if(!is.null(restrict))
  {
    if(restrict == "numeric")
    {
      isNumeric <- sapply(dataTable, is.numeric)
      dataTable <- dataTable[, isNumeric, drop = FALSE]
    } else if(restrict == "integer")
    {
      isInteger <- sapply(dataTable, is.integer)
      dataTable <- dataTable[, isInteger, drop = FALSE]
    }
  }

  # Only return independent variables in the table.
  dropColumns <- na.omit(match(c("primary", "class"), colnames(dataTable)))
  if(length(dropColumns) > 0) dataTable <- dataTable[, -dropColumns]
  
  if(!is.null(classes))
    list(dataTable = dataTable, classes = classes)
  else
    dataTable
}

.checkVariablesAndSame <- function(trainingMatrix, testingMatrix)
{
  if(ncol(trainingMatrix) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else if(ncol(trainingMatrix) != ncol(testingMatrix))
    stop("Training data set and testing data set contain differing numbers of features.")  
}

.doSelection <- function(measurements, classes, featureSets, metaFeatures, training, selectParams, trainParams,
                         predictParams, verbose)
{
  initialClass <- class(measurements)
  if(!"list" %in% initialClass)
    measurements <- list(data = measurements)  

  names(classes) <- rownames(measurements[[1]]) # In case training specified by sample IDs rather than numeric indices.
  trainClasses <- droplevels(classes[training])
  rankedSelected <- lapply(measurements, function(measurementsVariety)
  {
    if(is.function(selectParams@featureSelection))
    {
      paramList <- list(measurementsVariety[training, , drop = FALSE], trainClasses, verbose = verbose)
      selectFormals <- names(.methodFormals(selectParams@featureSelection))
      if("trainParams" %in% selectFormals) # Needs training and prediction functions for resubstitution error rate calculation.
        paramList <- append(paramList, c(trainParams = trainParams, predictParams = predictParams))
      if(!is.null(metaFeatures))
        paramList <- append(paramList, c(metaFeatures = metaFeatures[training, , drop = FALSE]))
      if("featureSets" %in% selectFormals) # Pass the sets on from runTest.
        paramList <- append(paramList, c(featureSets = featureSets))
      paramList <- append(paramList, c(selectParams@otherParams, datasetName = "N/A", selectionName = "N/A"))
      selection <- do.call(selectParams@featureSelection, paramList)

      if(class(selection) == "SelectResult")
      {
        if(length(selection@rankedFeatures) == 0) ranked <- numeric() else ranked <- selection@rankedFeatures[[1]]
        list(ranked, selection@chosenFeatures[[1]])
      } else { # List of such results for varieties.
        list(lapply(selection, function(variety) {
                               if(length(variety@rankedFeatures) == 0)
                                 numeric()
                               else
                                 variety@rankedFeatures[[1]]
                               }),
             lapply(selection, function(variety) variety@chosenFeatures[[1]]))
      }
    } else { # It is a list of functions for ensemble selection.
      featuresLists <- mapply(function(selector, selParams)
      {
        paramList <- list(measurementsVariety[training, , drop = FALSE], trainClasses, trainParams = trainParams,
                          predictParams = predictParams, verbose = verbose)
        paramList <- append(paramList, c(selParams, datasetName = "N/A", selectionName = "N/A"))
        do.call(selector, paramList)
      }, selectParams@featureSelection, selectParams@otherParams, SIMPLIFY = FALSE)

      if(class(featuresLists[[1]]) == "SelectResult") # No varieties were returned by the classifier used for resubstitution.
      {
        if(is.vector(featuresLists[[1]]@chosenFeatures[[1]])) # Data set is not MultiAssayExperiment, only variable ID tracked.
        {
          featuresCounts <- table(unlist(lapply(featuresLists, function(featureSet) featureSet@chosenFeatures[[1]])))
          selectedFeatures <- names(featuresCounts)[featuresCounts >= selectParams@minPresence]
          selectedFeatures
        } else { # Selected feature information is stored in a data frame.
          chosenFeaturesEnsemble <- do.call(rbind, lapply(featuresLists, function(featureSet) featureSet@chosenFeatures[[1]]))
          selectedFeatures <- chosenFeaturesEnsemble[plyr::count(chosenFeaturesEnsemble)[, "freq"] >= selectParams@minPresence, c("dataset", "feature")]
          selectedFeatures
        }
      } else { # The prediction function used for resubstitution returned a variety of lists.
        selectedFeatures <- lapply(1:length(featuresLists[[1]]), function(variety)
        {
          varietyFeatures <- lapply(featuresLists, function(selectList) selectList[[variety]]@chosenFeatures[[1]])
          if(is.vector(featuresLists[[1]]@chosenFeatures[[1]])) # Data set is not MultiAssayExperiment, only variable ID tracked.
          {
            featuresCounts <- table(unlist(varietyFeatures))
            selectedFeatures <- names(featuresCounts)[featuresCounts >= selectParams@minPresence]
            selectedFeatures
          } else { # Selected feature information is stored in a data frame.
            chosenFeaturesEnsemble <- do.call(rbind, varietyFeatures)
            selectedFeatures <- chosenFeaturesEnsemble[plyr::count(chosenFeaturesEnsemble)[, "freq"] >= selectParams@minPresence, c("dataset", "feature")]
            selectedFeatures
          }
        })
        names(selectedFeatures) <- names(featuresLists[[1]]) # Add variety names to selected list.
      }
      list(NULL, selectedFeatures)
    }
  })

  if(!"list" %in% initialClass) rankedSelected <- rankedSelected[[1]]
  rankedSelected
}

.doTransform <- function(measurements, training, transformParams, verbose)
{
  initialClass <- class(measurements)
  if(!"list" %in% class(measurements))
    measurements <- list(data = measurements)
  
  transformed <- lapply(measurements, function(measurementsVariety)
  {
    paramList <- list(measurementsVariety, training = training) # Often a central point, like the mean, is used for subtraction or standardisation of values. Pass this to the transformation function.
    if(length(transformParams@otherParams) > 0)
      paramList <- c(paramList, transformParams@otherParams)
    paramList <- c(paramList, verbose = verbose)
    do.call(transformParams@transform, paramList)
  })
  
  if(!"list" %in% initialClass) transformed <- transformed[[1]]
  if("list" %in% class(transformed[[1]])) transformed <- unlist(transformed, recursive = FALSE)
  transformed
}

.doTrain <- function(measurements, classes, training, testing, trainParams,
                            predictParams, verbose)
  # Re-use inside feature selection.
{
  initialClass <- class(measurements)
  if(!"list" %in% class(measurements)) # Will be a DataFrame.
    measurements <- list(data = measurements)

  names(classes) <- rownames(measurements[[1]]) # In case training or testing specified by sample IDs rather than numeric indices.
  trainClasses <- droplevels(classes[training])
  trained <- mapply(function(measurementsVariety, variety)
  {
    measurementsTrain <- measurementsVariety[training, , drop = FALSE]
    measurementsTest <- measurementsVariety[testing, , drop = FALSE]
    if(variety != "data") # Single measurements table is in a list with name 'data'.
    {
      multiplierParams <- sapply(strsplit(variety, ",")[[1]], strsplit, split = '=')
      individiualParams <- lapply(multiplierParams, '[', 2)
      names(individiualParams) <- sapply(multiplierParams, '[', 1)
      individiualParams <- lapply(individiualParams, function(param) tryCatch(as.numeric(param), warning = function(warn){param}))
      trainFormals <- names(.methodFormals(trainParams@classifier))
      predictFormals <- character()
      if(!is.null(predictParams@predictor)) # Only check for formals if a function was specified by the user.
        predictFormals <- names(.methodFormals(predictParams@predictor))
      changeTrain <- intersect(names(individiualParams), trainFormals)
      changePredict <- intersect(names(individiualParams), predictFormals)
      trainParams@otherParams[changeTrain] <- individiualParams[changeTrain]
      predictParams@otherParams[changePredict] <- individiualParams[changePredict]
    }
    
    tuneIndex <- which(names(trainParams@otherParams) == "tuneParams")
    if(length(tuneIndex) > 0)
    {
      tuneCombinations <- expand.grid(trainParams@otherParams[[tuneIndex]])
      trainParams@otherParams <- trainParams@otherParams[-tuneIndex]
    } else tuneCombinations <- NULL

    if(trainParams@classifier@generic != "previousTrained")
      paramList <- list(measurementsTrain, trainClasses)
    else # Don't pass the measurements and classes, because a pre-existing classifier is used.
      paramList <- list()
    if(!is.null(predictParams@predictor)) # Training and prediction are separate.
    {
      if(is.null(tuneCombinations))
      {
        if(length(trainParams@otherParams) > 0)
          paramList <- c(paramList, trainParams@otherParams)
        paramList <- c(paramList, verbose = verbose)
        trained <- do.call(trainParams@classifier, paramList)
        if(verbose >= 2)
          message("Training completed.")  
        returnResult <- trained
      } else {# Tuning Parameter selection.
        trainedList <- list()
        tuneOptimise <- trainParams@otherParams[["tuneOptimise"]]
        trainParams@otherParams <- trainParams@otherParams[-match("tuneOptimise", names(trainParams@otherParams))] # Don't pass the tuning optimisation parameters directly to the classifier.
        performances <- apply(tuneCombinations, 1, function(tuneCombination)
        {
          tuneParams <- as.list(tuneCombination)
          if(length(trainParams@otherParams) > 0)
            paramList <- c(paramList, trainParams@otherParams)
          paramList <- c(paramList, tuneParams, verbose = verbose)
          trained <- do.call(trainParams@classifier, paramList)
          initialTrainClass <- class(trained)
          if(! "list" %in% initialTrainClass) trained <- list(trained)

          if("tune" %in% names(attributes(trained[[1]])))
              tuneParams <- c(tuneParams, attr(trained[[1]], "tune"))
          if(variety != "data")
            names(trained) <- paste(variety, names(trained), sep = ',')
          trainedList <<- c(trainedList, trained)
    
          lapply(trained, function(model)
          {        
            paramList <- list(model, measurementsTrain) # Model and test set same as training set.
            if(length(predictParams@otherParams) > 0)
              paramList <- c(paramList, predictParams@otherParams)
            paramList <- c(paramList, verbose = verbose)
            predicted <- do.call(predictParams@predictor, paramList)
            
            if(class(predicted) != "list" || sum(grepl('=', names(predicted))) == 0)
              predicted <- list(predicted)
            
            if(class(predicted[[1]]) == "data.frame") # Predictor returned both scores and classes; just use classes.
              predicted <- lapply(predicted, function(variety) variety[, sapply(variety, class) == "factor"])
            if(is.numeric(class(predicted[[1]]))) # Can't automatically decide on a threshold. Stop processing.
               stop("Only numeric predictions are available. Predicted classes must be provided.")

            lapply(predicted, function(predictions)
            {
              calcExternalPerformance(trainClasses, predictions, tuneOptimise[1])
            })
          })
        })
        
        chosenModels <- lapply(1:length(performances[[1]]), function(trainVariety)
        {
          lapply(1:length(trainVariety[[1]]), function(predictVariety)
          {
              performanceValues <- sapply(performances, function(tuneLevel) tuneLevel[[trainVariety]][[predictVariety]])
              if(tuneOptimise[2] == "lower")
                chosenTune <- which.min(performanceValues)[1]
              else # It is "higher"
                chosenTune <- which.max(performanceValues)[1]

              chosenModel <- trainedList[[chosenTune]]
              tuneParameters <- as.list(tuneCombinations[chosenTune, , drop = FALSE])
              # Concatenate in case the classifier interally did other parameter tuning and recorded it in the tune attribute.
              attr(chosenModel, "tune") <- c(attr(chosenModel, "tune"), tuneParameters)
              chosenModel
          })
        })
        
        chosenNames <- paste(rep(names(performances[[1]]), each = length(performances[[1]][[1]])), rep(names(performances[[1]][[1]]), length(performances[[1]])), sep = ',')
        chosenNames <- gsub("^,|,$", '', chosenNames)
        if(variety != "data")
          chosenNames <- paste(variety, chosenNames, sep = ',')
        if(length(chosenNames) == 0)
        {
          chosenModels <- chosenModels[[1]]
          if(length(chosenModels) == 1)
            chosenModels <- chosenModels[[1]]
        } else {
          chosenModels <- unlist(chosenModels, recursive = FALSE)
          names(chosenModels) <- chosenNames
        }
        if(verbose >= 2)
          message("Parameter tuning and training completed.")
  
        returnResult <- chosenModels
      }
    } else { # Some classifiers do training and testing with a single function.
      paramList <- append(paramList, list(measurementsTest))
      if(length(trainParams@otherParams) > 0)
      {
        whichList <- which(sapply(trainParams@otherParams, is.list))
        extras <- trainParams@otherParams
        if(length(whichList) > 0) # Used when selected features is passed as an intermediate and there are multiple varieties.
        {
          extras[whichList] <- lapply(whichList, function(paramIndex)
                               {
                                 varietyIndex <- match(variety, names(extras[[paramIndex]]))
                                 if(!is.na(varietyIndex))
                                   extras[[paramIndex]][[varietyIndex]]
                                 else
                                   extras[[whichList]]
                               })
        }
        paramList <- c(paramList, extras)
      }
      if(length(predictParams@otherParams) > 0)
      {
        whichList <- which(sapply(predictParams@otherParams, is.list))
        extras <- predictParams@otherParams
        if(length(whichList) > 0) # Used when selected features is passed as an intermediate and there are multiple varieties.
        {
          extras[whichList] <- lapply(whichList, function(paramIndex)
                               {
                                 varietyIndex <- match(variety, names(extras[[paramIndex]]))
                                 if(!is.na(varietyIndex))
                                   extras[[paramIndex]][[varietyIndex]]
                                 else
                                   extras[[whichList]]
                               })
        }
        paramList <- c(paramList, extras)
      }
      paramList <- c(paramList, verbose = verbose)
            
      if(is.null(tuneCombinations))
      {
        predictions <- do.call(trainParams@classifier, paramList)
        if(verbose >= 2)
          message("Training and prediction completed.")    
        returnResult <- predictions
      } else { # Tuning Parameter selection.
        tuneOptimise <- trainParams@otherParams[["tuneOptimise"]]
        trainParams@otherParams <- trainParams@otherParams[-match("tuneOptimise", names(trainParams@otherParams))] # Don't pass the tuning optimisation parameters directly to the classifier.
        performances <- apply(tuneCombinations, 1, function(tuneCombination)
        {
          tuneParams <- as.list(tuneCombination)
          names(tuneParams) <- colnames(tuneCombinations)
          paramList <- c(paramList, tuneParams)
          trained <- do.call(trainParams@classifier, paramList)
          if(class(trained) != "list" || sum(grepl('=', names(trained))) == 0)
            predictedClasses <- list(trained)
          else
            predictedClasses <- trained
          if(class(trained[[1]]) == "data.frame") # Predictor returned both scores and classes. Just use classes.
            predictedFactor <- lapply(trained, function(variety) variety[, sapply(variety, class) == "factor"])
          else if(class(trained[[1]]) == "factor")
            predictedFactor <- trained
          else if(is.numeric(class(trained[[1]]))) # Can't automatically decide on a threshold. Stop processing.
            stop("Only numeric predictions are available. Predicted classes must be provided.")
  
          tunePredictions <- lapply(1:length(predictedFactor), function(predictIndex)
          {
            performanceValue <- calcExternalPerformance(trainClasses, predictedFactor[[predictIndex]],
                                                        tuneOptimise[1])
            list(predictedClasses[[predictIndex]], performanceValue)
          })
        })

        chosenPredictions <- lapply(1:length(performances[[1]]), function(predictVariety)
        {
          performanceValues <- sapply(performances, function(tuneLevel) tuneLevel[[predictVariety]][[2]]) # Value is in second position.
          if(tuneOptimise[2] == "lower")
            chosenTune <- which.min(performanceValues)[1]
          else # It is "higher"
            chosenTune <- which.max(performanceValues)[1]
            
          chosenPredict <- performances[[chosenTune]][[predictVariety]][[1]] # Prediction object is in position 1.
          attr(chosenPredict, "tune") <- as.list(tuneCombinations[chosenTune, , drop = FALSE])
        })
        names(chosenPredictions) <- names(trained)
        if(variety != "data")
          names(chosenPredictions) <- paste(variety, names(trained), sep = '')
        if(length(chosenPredictions) == 1) chosenPredictions <- chosenPredictions[[1]]
        
        if(verbose >= 2)
          message("Parameter tuning and classification completed.")
        returnResult <- chosenPredictions
      }
    }
    returnResult
  }, measurements, names(measurements), SIMPLIFY = FALSE)

  if(!"list" %in% initialClass) trained <- trained[[1]]
  if(!isS4(trained) && "list" %in% class(trained[[1]])) 
  {
    trainNames <- sapply(trained, names)
    varietyNames <- sapply(trained[[1]], names)
    trained <- unlist(trained, recursive = FALSE)
    names(trained) <- paste(rep(trainNames, each = length(varietyNames), rep(varietyNames, length(trainNames))), sep = ',')
  }
  
  trained
}

.doTest <- function(trained, measurements, testing, predictParams, verbose)
  # Re-use inside feature selection.
{
  initialClass <- class(measurements)
  if(!"list" %in% initialClass)
  {
    trained <- list(model = trained)
    measurements <- list(data = measurements)
  }

  predicted <- mapply(function(model, data, variety)
  {
    if(!is.null(predictParams@predictor))
    {
      testMeasurements <- data[testing, , drop = FALSE]
      if(variety != "data") # Single measurements table is in a list with name 'data'.
      {
        multiplierParams <- sapply(strsplit(variety, ",")[[1]], strsplit, split = '=')
        individiualParams <- lapply(multiplierParams, '[', 2)
        names(individiualParams) <- sapply(multiplierParams, '[', 1)
        individiualParams <- lapply(individiualParams, function(param) tryCatch(as.numeric(param), warning = function(warn){param}))
        change <- intersect(names(individiualParams), names(predictParams@otherParams))
        predictParams@otherParams[change] <- individiualParams[change]
      }      
    
      paramList <- list(model, testMeasurements)
      if(length(predictParams@otherParams) > 0)
        paramList <- c(paramList, predictParams@otherParams)
      paramList <- c(paramList, verbose = verbose)
      prediction <- do.call(predictParams@predictor, paramList)
    } else
    {
      prediction <- model
    }
    
    if(verbose >= 2)
      message("Prediction completed.")    
    prediction
  }, trained, measurements, names(measurements), SIMPLIFY = FALSE)
  
  if(!"list" %in% initialClass) predicted <- predicted[[1]]
  if("list" %in% class(predicted[[1]])) predicted <- unlist(predicted, recursive = FALSE)
  predicted
}

.pickFeatures <- function(measurements, classes, featureSets, datasetName, trainParams, predictParams,
                          resubstituteParams, ordering, selectionName, verbose)
{
  maxFeatures <- max(resubstituteParams@nFeatures)
  if(maxFeatures > ncol(measurements))
    stop("Feature selection specified to consider as many as ", maxFeatures, " features, but data set has only ", ncol(measurements), " features.")

  if(is.null(featureSets))
  {
    orderedList <- as.list(ordering[1:maxFeatures])
  } else { # Group by sets.
    if(!"Pairs" %in% class(featureSets))
      orderedList <- split(1:ncol(measurements), S4Vectors::mcols(measurements)[["original"]])[ordering[1:maxFeatures]]
    else # Is a Pairs object.
      orderedList <- featureSets[ordering[1:maxFeatures]]
  }

  performances <- sapply(resubstituteParams@nFeatures, function(topFeatures)
  {
    if(!"Pairs" %in% class(featureSets))
    {
      measurementsSubset <- measurements[, unlist(orderedList[1:topFeatures]), drop = FALSE]
      trained <- .doTrain(measurementsSubset, classes, 1:nrow(measurementsSubset), 1:nrow(measurementsSubset),
                          trainParams, predictParams, verbose)
    } else { # Pairs; don't subset.
      trainParams@otherParams <- c(trainParams@otherParams, featurePairs = featureSets[1:topFeatures])
      trained <- .doTrain(measurements, classes, 1:nrow(measurements), 1:nrow(measurements),
                          trainParams, predictParams, verbose)
    }

    if(!is.null(predictParams@predictor))
    {
        predictions <- .doTest(trained, measurementsSubset, 1:nrow(measurementsSubset),
                               predictParams, verbose)
    } else {
      predictions <- trained
    }
    
    if("list" %in% class(predictions)) # Mutiple varieties of predictions.
    {
      if(class(predictions[[1]]) == "data.frame")
        predictedClasses <- lapply(predictions, function(set) set[, sapply(set, class) == "factor"])
      else
        predictedClasses <- predictions
    } else { # A single variety of prediction.
      if(class(predictions) == "data.frame")
        predictedClasses <- predictions[, sapply(predictions, class) == "factor"]
      else
        predictedClasses <- predictions
    }

    if("list" %in% class(predictedClasses))
    {
      lapply(predictedClasses, function(classSet) calcExternalPerformance(classes, classSet, resubstituteParams@performanceType))
    } else {
      calcExternalPerformance(classes, predictedClasses, resubstituteParams@performanceType)
    }
  })
  
  if(class(performances) == "numeric")
    performances <- matrix(performances, ncol = length(performances), byrow = TRUE)

  pickedFeatures <- apply(performances, 1, function(varietyPerformances)
                    {
                      if(resubstituteParams@better == "lower")
                        1:(resubstituteParams@nFeatures[which.min(varietyPerformances)[1]])     
                      else
                        1:(resubstituteParams@nFeatures[which.max(varietyPerformances)[1]])
                    })

  if(is.matrix(pickedFeatures)) # Same number of features picked for each variety. Coerce to list.
    pickedFeatures <- as.list(as.data.frame(pickedFeatures))
  
  if(verbose == 3)
    message("Features selected.")

  rankedFeatures <- lapply(1:length(pickedFeatures), function(variety) ordering)
  pickedFeatures <- lapply(pickedFeatures, function(pickedSet) ordering[pickedSet])
  
  if(!is.null(S4Vectors::mcols(measurements)) && "dataset" %in% colnames(S4Vectors::mcols(measurements))) # Table describing source table and variable name is present.
  {
    varInfo <- S4Vectors::mcols(measurements)
    rankedFeatures <- lapply(rankedFeatures, function(features) varInfo[features, ])
    pickedFeatures <- lapply(pickedFeatures, function(features) varInfo[features, ])
  } else { # Vectors of feature names.
    if(is.null(featureSets))
    {
      rankedFeatures <- lapply(rankedFeatures, function(features) colnames(measurements)[features])
      pickedFeatures <- lapply(pickedFeatures, function(features) colnames(measurements)[features])
    } else {
      if(!"Pairs" %in% class(featureSets))
      {
        rankedFeatures <- lapply(rankedFeatures, function(features) names(featureSets@sets)[features])
        pickedFeatures <- lapply(pickedFeatures, function(features) names(featureSets@sets)[features])
      } else {
        rankedFeatures <- lapply(rankedFeatures, function(features) featureSets[features])
        pickedFeatures <- lapply(pickedFeatures, function(features) featureSets[features])
      }
    }
  }

  if(is.null(featureSets))
    totalFeatures <- ncol(measurements)
  else if(!"Pairs" %in% class(featureSets))
    totalFeatures <- length(featureSets@sets)
  else
    totalFeatures <- length(featureSets)
  
  selectResults <- lapply(1:length(rankedFeatures), function(variety)
  {
    SelectResult(datasetName, selectionName, totalFeatures,
                 list(rankedFeatures[[variety]]), list(pickedFeatures[[variety]]))
  })
  names(selectResults) <- names(pickedFeatures)

  if(length(selectResults) == 1) selectResults <- selectResults[[1]] else selectResults
}

.validationText <- function(result)
{
  switch(result@validation[[1]],
  permuteFold = paste(result@validation[[2]], "Permutations,", result@validation[[3]], "Folds"),
  fold = paste(result@validation[[2]], "-fold cross-validation", sep = ''),
  leave = paste("Leave", result@validation[[2]], "Out"),
  split = paste(result@validation[[2]], "Permutations,", result@validation[[3]], "% Test"),
  independent = "Independent Set")
}

.binValues <- function(values, nBins)
{
  ordering <- order(values)
  binID <- rep(1:nBins, each = length(values) / nBins)
  if(length(binID) < length(values))
    binID <- c(binID, rep(max(binID) + 1, length(values) - length(binID)))
  bins <- split(ordering, binID)  
  binID <- numeric()
  binID[unlist(bins)] <- rep(as.numeric(names(bins)), sapply(bins, length))
  binID
}

.methodFormals <- function(f) {
  tryCatch({
    fdef <- getGeneric(f)
    method <- selectMethod(fdef, "DataFrame")
    genFormals <- base::formals(fdef)
    b <- body(method)
    if(is(b, "{") && is(b[[2]], "<-") && identical(b[[2]][[2]], as.name(".local"))) {
      local <- eval(b[[2]][[3]])
      if(is.function(local))
        return(formals(local))
      warning("Expected a .local assignment to be a function. Corrupted method?")
    }
    genFormals
  },
    error = function(error) {
      formals(f)
    })
}

.densitiesCrossover <- function(densities) # A list of densities created by splinefun.
{
  
  if(!all(table(unlist(lapply(densities, function(density) density[['x']]))) == length(densities)))
    stop("x positions are not the same for all of the densities.")
  
  lapply(1:length(densities), function(densityIndex) # All crossing points with other class densities.
  {
    unlist(lapply(setdiff(1:length(densities), densityIndex), function(otherIndex)
    {
      allDifferences <- densities[[densityIndex]][['y']] - densities[[otherIndex]][['y']]
      crosses <- which(diff(sign(allDifferences)) != 0)
      crosses <- sapply(crosses, function(cross) # Refine location for plateaus.
      {
        isSmall <- rle(allDifferences[(cross+1):length(allDifferences)] < 0.000001)
        if(isSmall[["values"]][1] == "TRUE")
          cross <- cross + isSmall[["lengths"]][1] / 2
        cross
      })
      if(length(crosses) > 1 && densities[[densityIndex]][['y']][crosses[1]] < 0.000001 && densities[[densityIndex]][['y']][crosses[length(crosses)]] < 0.000001)
        crosses <- crosses[-c(1, length(crosses))] # Remove crossings at ends of densities.      
      densities[[densityIndex]][['x']][crosses]
    }))
  })
}

.dlda <- function(x, y, prior = NULL){ # Remove once sparsediscrim is reinstated to CRAN.
  obj <- list()
  obj$labels <- y
  obj$N <- length(y)
  obj$p <- ncol(x)
  obj$groups <- levels(y)
  obj$num_groups <- nlevels(y)

  est_mean <- "mle"

  # Error Checking
  if (!is.null(prior)) {
    if (length(prior) != obj$num_groups) {
      stop("The number of 'prior' probabilities must match the number of classes in 'y'.")
    }
    if (any(prior <= 0)) {
      stop("The 'prior' probabilities must be nonnegative.")
    }
    if (sum(prior) != 1) {
      stop("The 'prior' probabilities must sum to one.")
    }
  }
  if (any(table(y) < 2)) {
    stop("There must be at least 2 observations in each class.")
  }

  # By default, we estimate the 'a priori' probabilities of class membership with
  # the MLEs (the sample proportions).
  if (is.null(prior)) {
    prior <- as.vector(table(y) / length(y))
  }

  # For each class, we calculate the MLEs (or specified alternative estimators)
  # for each parameter used in the DLDA classifier. The 'est' list contains the
  # estimators for each class.
  obj$est <- tapply(seq_along(y), y, function(i) {
    stats <- list()
    stats$n <- length(i)
    stats$xbar <- colMeans(x[i, , drop = FALSE])
    stats$var <- with(stats, (n - 1) / n * apply(x[i, , drop = FALSE], 2, var))
    stats
  })

  # Calculates the pooled variance across all classes.
  obj$var_pool <- Reduce('+', lapply(obj$est, function(x) x$n * x$var)) / obj$N

  # Add each element in 'prior' to the corresponding obj$est$prior
  for(k in seq_len(obj$num_groups)) {
    obj$est[[k]]$prior <- prior[k]
  }
  class(obj) <- "dlda"
  obj
}

.predict <- function(object, newdata, ...) { # Remove once sparsediscrim is reinstated to CRAN.
  if (!inherits(object, "dlda"))  {
    stop("object not of class 'dlda'")
  }
  if (is.vector(newdata)) {
    newdata <- as.matrix(newdata)
  }

  scores <- apply(newdata, 1, function(obs) {
    sapply(object$est, function(class_est) {
      with(class_est, sum((obs - xbar)^2 / object$var_pool) + log(prior))
    })
  })

  if (is.vector(scores)) {
    min_scores <- which.min(scores)
  } else {
    min_scores <- apply(scores, 2, which.min)
  }

  # Posterior probabilities via Bayes Theorem
  means <- lapply(object$est, "[[", "xbar")
  covs <- replicate(n=object$num_groups, object$var_pool, simplify=FALSE)
  priors <- lapply(object$est, "[[", "prior")
  posterior <- .posterior_probs(x=newdata,
                               means=means,
                               covs=covs,
                               priors=priors)

  class <- factor(object$groups[min_scores], levels = object$groups)

  list(class = class, scores = scores, posterior = posterior)
}

.posterior_probs <- function(x, means, covs, priors) { # Remove once sparsediscrim is reinstated to CRAN.
  if (is.vector(x)) {
    x <- matrix(x, nrow=1)
  }
  x <- as.matrix(x)

  posterior <- mapply(function(xbar_k, cov_k, prior_k) {
    if (is.vector(cov_k)) {
      post_k <- apply(x, 1, function(obs) {
        .dmvnorm_diag(x=obs, mean=xbar_k, sigma=cov_k)
      })
    } else {
      post_k <- dmvnorm(x=x, mean=xbar_k, sigma=cov_k)
    }
    prior_k * post_k
  }, means, covs, priors)

  if (is.vector(posterior)) {
    posterior <- posterior / sum(posterior)
  } else {
    posterior <- posterior / rowSums(posterior)
  }

  posterior
}

.dmvnorm_diag <- function(x, mean, sigma) { # Remove once sparsediscrim is reinstated to CRAN.
  exp(sum(dnorm(x, mean=mean, sd=sqrt(sigma), log=TRUE)))
}
