setOldClass("pamrtrained")

.doSelection <- function(expression, training, selectParams, trainParams,
                         predictParams, verbose)
{
  initialClass <- class(expression)
  if(class(expression) != "list")
    expression <- list(data = expression)  

  rankedSelected <- lapply(expression, function(expressionVariety)
  {
    if(is.function(selectParams@featureSelection))
    {
      paramList <- list(expressionVariety[, training], verbose = verbose)
      if(selectParams@featureSelection@generic[1] != "previousSelection")
        paramList <- append(paramList, c(trainParams = trainParams, predictParams = predictParams))
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
        paramList <- list(expressionVariety[, training], trainParams = trainParams,
                          predictParams = predictParams, verbose = verbose)
        paramList <- append(paramList, selParams)
        do.call(selector, paramList)
      }, selectParams@featureSelection, selectParams@otherParams, SIMPLIFY = FALSE)
      
      if(class(featuresLists[[1]]) == "SelectResult")
      {
        featuresCounts <- table(unlist(lapply(featuresLists, function(featureSet) featureSet@chosenFeatues[[1]])))
        selectedFeatures <- as.integer(names(featuresCounts))[featuresCounts >= selectParams@minPresence]
      } else { # The prediction function used for resubstitution returned a variety of lists.
        selectedFeatures <- lapply(1:length(featuresLists[[1]]), function(variety)
        {
          varietyFeatures <- lapply(featuresLists, function(selectList) selectList[[variety]]@chosenFeatues[[1]])
          featuresCounts <- table(unlist(varietyFeatures))
          as.integer(names(featuresCounts))[featuresCounts >= selectParams@minPresence]
        })
      }
      list(NULL, selectedFeatures)
    }
  })

  if(initialClass != "list") rankedSelected <- rankedSelected[[1]]
  rankedSelected
}

.doTransform <- function(expression, transformParams, verbose)
{
  initialClass <- class(expression)
  if(class(expression) != "list")
    expression <- list(data = expression)
  
  transformed <- lapply(expression, function(expressionVariety)
  {
    paramList <- list(expressionVariety)
    if(length(transformParams@otherParams) > 0)
      paramList <- c(paramList, transformParams@otherParams)
    paramList <- c(paramList, verbose = verbose)
    do.call(transformParams@transform, paramList)
  })
  
  if(initialClass != "list") transformed <- transformed[[1]]
  if(class(transformed[[1]]) == "list") transformed <- unlist(transformed, recursive = FALSE)
  transformed
}

.doTrain <- function(expression, training, testing, trainParams,
                            predictParams, verbose)
  # Re-use inside feature selection.
{
  initialClass <- class(expression)
  if(class(expression) != "list")
    expression <- list(data = expression)
  
  trained <- mapply(function(expressionVariety, variety)
  {
    classes <- pData(expressionVariety)[, "class"]
    expressionTrain <- exprs(expressionVariety[, training])
    expressionTest <- exprs(expressionVariety[, testing])
    if(trainParams@transposeExpression == TRUE) # For classifiers that consider columns as samples and genes as rows.
    {
      expressionTrain <- t(expressionTrain)
      expressionTest <- t(expressionTest)
    }
    
    if(variety != "data") # Single expression set is in a list with name 'data'.
    {
      multiplierParams <- sapply(strsplit(variety, ",")[[1]], strsplit, split = '=')
      individiualParams <- lapply(multiplierParams, '[', 2)
      names(individiualParams) <- sapply(multiplierParams, '[', 1)
      individiualParams <- lapply(individiualParams, function(param) tryCatch(as.numeric(param), warning = function(warn){param}))
      changeTrain <- intersect(names(individiualParams), names(trainParams@otherParams))
      changePredict <- intersect(names(individiualParams), names(predictParams@otherParams))
      trainParams@otherParams[changeTrain] <- individiualParams[changeTrain]
      predictParams@otherParams[changePredict] <- individiualParams[changePredict]
    }
    
    tuneIndex <- which(names(trainParams@otherParams) == "tuneParams")
    if(length(tuneIndex) > 0)
    {
      tuneCombinations <- expand.grid(trainParams@otherParams[[tuneIndex]])
      trainParams@otherParams <- trainParams@otherParams[-tuneIndex]
    } else tuneCombinations <- NULL
    paramList <- list(expressionTrain, classes[training])
    if(trainParams@doesTests == FALSE) # Training and prediction are separate.
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
      } else { # Tuning Parameter selection.
        trainedList <- list()
        resubstituteParams <- trainParams@otherParams[["resubstituteParams"]]
        performances <- apply(tuneCombinations, 1, function(tuneCombination)
        {
          tuneParams <- as.list(tuneCombination)
          names(tuneParams) <- colnames(tuneCombinations)
          if(length(trainParams@otherParams) > 0)
            paramList <- c(paramList, trainParams@otherParams)
          paramList <- c(paramList, tuneParams, verbose = verbose)
          trained <- do.call(trainParams@classifier, paramList)
          initialTrainClass <- class(trained)
          if(initialTrainClass != "list") trained <- list(trained)
          if(variety != "data")
            names(trained) <- paste(variety, names(trained), sep = ',')
          trainedList <<- c(trainedList, trained)
          
          lapply(trained, function(model)
          {            
            paramList <- list(model, expressionTrain) # Model and test set same as training set.
            if(length(predictParams@otherParams) > 0)
              paramList <- c(paramList, predictParams@otherParams)
            paramList <- c(paramList, tuneParams, verbose = verbose)
            predicted <- do.call(predictParams@predictor, paramList)
            if(class(predicted) != "list" || sum(grepl('=', names(predicted))) == 0)
              predicted <- list(predictParams@getClasses(predicted))
            else
              predicted <- lapply(predicted, predictParams@getClasses)
            if(class(predicted[[1]]) == "data.frame") # Predictor returned both scores and labels. Just use labels.
              predicted <- lapply(predicted, function(variety) variety[, sapply(variety, class) == "factor"])
            if(is.numeric(class(predicted[[1]]))) # Can't automatically decide on a threshold. Stop processing.
               stop("Only numeric predictions are available. Predicted class labels must be provided.")
            lapply(predicted, function(predictions)
            {
              classData <- ROCR::prediction(as.numeric(predictions), factor(classes[training], ordered = TRUE))
              if(resubstituteParams@performanceType == "balanced")
              {
                  falseNegativeRate <- ROCR::performance(classData, "fnr")@y.values[[1]][2]
                  falsePositiveRate <- ROCR::performance(classData, "fpr")@y.values[[1]][2]
                  performanceValue <- mean(falseNegativeRate, falsePositiveRate)
              } else
              {
                performanceList <- append(list(classData, resubstituteParams@performanceType), resubstituteParams@otherParams)
                performanceData <- do.call(ROCR::performance, performanceList)
                performanceValue <- performanceData@y.values[[1]][2]
              }
                performanceValue
            })
          })
        })
        chosenModels <- lapply(1:length(performances[[1]]), function(trainVariety)
        {
          lapply(1:length(trainVariety[[1]]), function(predictVariety)
          {
              performanceValues <- sapply(performances, function(tuneLevel) tuneLevel[[trainVariety]][[predictVariety]])
              if(resubstituteParams@better == "lower")
                chosenTune <- which.min(performanceValues)
              else
                chosenTune <- which.max(performanceValues)

              chosenModel <- trainedList[[chosenTune]]
              tuneParameters <- as.list(tuneCombinations[chosenTune, , drop = FALSE])
              attr(tuneParameters, "out.attrs") <- NULL
              attr(chosenModel, "tune") <- tuneParameters
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
      paramList <- append(paramList, list(expressionTest))
      if(length(predictParams@otherParams) > 0)
        paramList <- c(paramList, predictParams@otherParams)
      paramList <- c(paramList, verbose = verbose)
            
      if(is.null(tuneCombinations))
      {
        predictions <- do.call(trainParams@classifier, paramList)
        if(verbose >= 2)
          message("Training and classification completed.")    
        returnResult <- predictions
      } else { # Tuning Parameter selection.
        resubstituteParams <- trainParams@otherParams[["resubstituteParams"]]
        performances <- apply(tuneCombinations, 1, function(tuneCombination)
        {
          tuneParams <- as.list(tuneCombination)
          names(tuneParams) <- colnames(tuneCombinations)
          paramList <- c(paramList, tuneParams)
          trained <- do.call(trainParams@classifier, paramList)
          if(class(trained) != "list" || sum(grepl('=', names(trained))) == 0)
            predictedClasses <- list(predictParams@getClasses(trained))
          else
            predictedClasses <- lapply(trained, predictParams@getClasses)
          if(class(trained[[1]]) == "data.frame") # Predictor returned both scores and labels. Just use labels.
            predictedFactor <- lapply(trained, function(variety) variety[, sapply(variety, class) == "factor"])
          else if(class(trained[[1]]) == "factor")
            predictedFactor <- trained
          else if(is.numeric(class(trained[[1]]))) # Can't automatically decide on a threshold. Stop processing.
          stop("Only numeric predictions are available. Predicted class labels must be provided.")
          
          tunePredictions <- lapply(1:length(predictedFactor), function(predictIndex)
          {
            classData <- ROCR::prediction(as.numeric(predictedFactor[[predictIndex]]), factor(classes[training], ordered = TRUE))
            
            if(resubstituteParams@performanceType == "balanced")
            {
              falseNegativeRate <- ROCR::performance(classData, "fnr")@y.values[[1]][2]
              falsePositiveRate <- ROCR::performance(classData, "fpr")@y.values[[1]][2]
              performanceValue <- mean(falseNegativeRate, falsePositiveRate)
            } else
            {
              performanceList <- append(list(classData, resubstituteParams@performanceType), resubstituteParams@otherParams)
              performanceData <- do.call(ROCR::performance, performanceList)
              performanceValue <- performanceData@y.values[[1]][2]
            }
            list(predictedClasses[[predictIndex]], performanceValue)
          })
        })
        
        chosenPredictions <- lapply(1:length(performances[[1]]), function(predictVariety)
        {
          performanceValues <- sapply(performances, function(tuneLevel) tuneLevel[[predictVariety]][[2]]) # Value is in second position.
          if(resubstituteParams@better == "lower")
            chosenTune <- which.min(performanceValues)
          else
            chosenTune <- which.max(performanceValues)
            
          chosenPredict <- performances[[chosenTune]][[predictVariety]][[1]] # Prediction object is in position 1.
          attr(chosenPredict, "tune") <- as.list(tuneCombinations[chosenTune, , drop = FALSE])
        })
        names(chosenPredictions) <- names(trained)
        if(variety != "data")
          names(chosenPredictions) <- paste(variety, names(trained), sep = '')
        if(length(chosenPredictions) == 1) chosenPredictions <- chosenPredictions[[1]]
        returnResult <- chosenPredictions
      }
    }
    returnResult
  }, expression, names(expression), SIMPLIFY = FALSE)

  if(initialClass != "list") trained <- trained[[1]]
  if(class(trained[[1]]) == "list") 
  {
    trainNames <- sapply(trained, names)
    varietyNames <- sapply(trained[[1]], names)
    trained <- unlist(trained, recursive = FALSE)
    names(trained) <- paste(rep(trainNames, each = length(varietyNames), rep(varietyNames, length(trainNames))), sep = ',')
  }
  
  trained
}

.doTest <- function(trained, expression, testing, predictParams, verbose)
  # Re-use inside feature selection.
{
  initialClass <- class(expression)
  if(initialClass != "list")
  {
    trained <- list(model = trained)
    expression <- list(data = expression)
  }
  
  predicted <- mapply(function(model, data, variety)
  {
    if(!grepl("{}", paste(capture.output(predictParams@predictor), collapse = ''), fixed = TRUE))
    {
      testExpression <- exprs(data)[, testing]
      if(predictParams@transposeExpression == TRUE)
        testExpression <- t(testExpression)
      
      if(variety != "data") # Single expression set is in a list with name 'data'.
      {
        multiplierParams <- sapply(strsplit(variety, ",")[[1]], strsplit, split = '=')
        individiualParams <- lapply(multiplierParams, '[', 2)
        names(individiualParams) <- sapply(multiplierParams, '[', 1)
        individiualParams <- lapply(individiualParams, function(param) tryCatch(as.numeric(param), warning = function(warn){param}))
        change <- intersect(names(individiualParams), names(predictParams@otherParams))
        predictParams@otherParams[change] <- individiualParams[change]
      }      
    
      paramList <- list(model, testExpression)
      if(length(predictParams@otherParams) > 0)
        paramList <- c(paramList, predictParams@otherParams)
      paramList <- c(paramList, verbose = verbose)
       prediction <- do.call(predictParams@predictor, paramList)
    } else
    {
      prediction <- model
    }
    
    if(class(prediction) != "list" || sum(grepl('=', names(prediction))) == 0)
      predictions <- predictParams@getClasses(prediction)
    else
      predictions <- lapply(prediction, predictParams@getClasses)
    
    if(verbose >= 2)
      message("Prediction completed.")    
    predictions
  }, trained, expression, names(expression), SIMPLIFY = FALSE)

  if(initialClass != "list") predicted <- predicted[[1]]
  if(class(predicted[[1]]) == "list") predicted <- unlist(predicted, recursive = FALSE)
  predicted
}

.pickRows <- function(expression, datasetName, trainParams, predictParams, resubstituteParams, orderedFeatures, selectionName, verbose)
{
  performances <- sapply(resubstituteParams@nFeatures, function(topFeatures)
  {
    expressionSubset <- expression[orderedFeatures[1:topFeatures], ]
    trained <- .doTrain(expressionSubset, 1:ncol(expressionSubset), 1:ncol(expressionSubset),
                        trainParams, predictParams, verbose)
    
    if(trainParams@doesTests == FALSE)
    {
        predictions <- .doTest(trained, expressionSubset, 1:ncol(expressionSubset),
                               predictParams, verbose)
    } else {
      if(class(trained) != "list" || sum(grepl('=', names(trained))) == 0)
        predictions <- predictParams@getClasses(trained)
      else
        predictions <- lapply(trained, function(model) .doTest(model, expressionSubset, 1:ncol(expressionSubset),
                                                               predictParams, verbose))
    }
    if(class(predictions) == "list") # Mutiple varieties of predictions.
    {
      if(class(predictions[[1]]) == "data.frame")
        labels <- lapply(predictions, function(set) set[, sapply(set, class) == "factor"])
      else
        labels <- predictions
      classData <- ROCR::prediction(lapply(labels, as.numeric), lapply(1:length(predictions), function(variety) factor(pData(expression)[, "class"], ordered = TRUE)))
    } else {
      if(class(predictions) == "data.frame")
        labels <- predictions[, sapply(predictions, class) == "factor"]
      else
        labels <- predictions
      classData <- ROCR::prediction(as.numeric(labels), factor(pData(expression)[, "class"], ordered = TRUE))
    }
    
    if(resubstituteParams@performanceType == "balanced")
    {
      falseNegativeRate <- sapply(ROCR::performance(classData, "fnr")@y.values, "[[", 2)
      falsePositiveRate <- sapply(ROCR::performance(classData, "fpr")@y.values, "[[", 2)
      performanceValues <- rowMeans(matrix(c(falseNegativeRate, falsePositiveRate), ncol = 2))
      if(length(performanceValues) > 1)
        names(performanceValues) <- names(predictions)
      performanceValues
    } else {
      performanceList <- append(list(classData, resubstituteParams@performanceType), resubstituteParams@otherParams)
      performanceData <- do.call(ROCR::performance, performanceList)
      performanceValues <- sapply(performanceData@y.values, "[[", 2)
      if(length(performanceValues) > 1)
        names(performanceValues) <- names(predictions)
      performanceValues
    }
  })
  
  if(class(performances) == "numeric")
    performances <- matrix(performances, ncol = length(performances), byrow = TRUE)

  pickedRows <- apply(performances, 1, function(varietyPerformances)
                {
                  if(resubstituteParams@better == "lower")
                    1:(resubstituteParams@nFeatures[which.min(varietyPerformances)[1]])     
                  else
                    1:(resubstituteParams@nFeatures[which.max(varietyPerformances)[1]])
                })

  if(is.matrix(pickedRows))
    pickedRows <- as.list(as.data.frame(pickedRows))
  
  if(verbose == 3)
    message("Features selected.")
  
  rankedFeatures <- lapply(1:length(pickedRows), function(variety) orderedFeatures)
  pickedFeatures <- lapply(pickedRows, function(pickedSet) orderedFeatures[pickedSet])
  
  selectResults <- lapply(1:length(rankedFeatures), function(variety)
  {
    SelectResult(datasetName, selectionName,
                 list(rankedFeatures[[variety]]), list(pickedFeatures[[variety]]))
  })
  names(selectResults) <- names(pickedRows)
  
  if(length(selectResults) == 1) selectResults <- selectResults[[1]] else selectResults
}

.validationText <- function(result)
{
  switch(result@validation[[1]],
  fold = paste(result@validation[[2]], "Resamples,", result@validation[[3]], "Folds"),
  leave = paste("Leave", result@validation[[2]], "Out"),
  split = paste(result@validation[[2]], "Resamples,", result@validation[[3]], "% Test"),
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

.methodFormals <- function(f, signature = character()) {
  fdef <- getGeneric(f)
  method <- selectMethod(fdef, signature)
  genFormals <- base::formals(fdef)
  b <- body(method)
  if(is(b, "{") && is(b[[2]], "<-") && identical(b[[2]][[2]], as.name(".local"))) {
    local <- eval(b[[2]][[3]])
    if(is.function(local))
      return(formals(local))
    warning("Expected a .local assignment to be a function. Corrupted method?")
  }
  genFormals
} 
