setOldClass("pamrtrained")

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
  }
  list(measurements = measurements, classes = classes)
}

.MAEtoWideTable <- function(measurements, targets, restrict = "numeric")
{
  if(!all(targets %in% c(names(measurements), "clinical")))
    stop("Some table names in 'targets' are not assay names in 'measurements' or \"clinical\".")  
  
  if("clinical" %in% targets)
  {
    clinicalColumns <- colnames(colData(measurements))
    targets <- targets[-match("clinical", targets)]
  } else {
    clinicalColumns <- NULL
  }
  measurements <- measurements[, , targets]
  dataTable <- wideFormat(measurements, colDataCols = clinicalColumns,
                                  check.names = FALSE)
  if("class" %in% colnames(dataTable))
    classes <- dataTable[, "class"]
  else
    classes <- NULL
  if(!is.null(restrict))
  {
    if(restrict == "numeric")
    {
      isNumeric <- apply(dataTable, 2, is.numeric)
      dataTable <- dataTable[, isNumeric, drop = FALSE]
    } else if(restrict == "integer")
    {
      isInteger <- apply(dataTable, 2, is.integer)
      dataTable <- dataTable[, isInteger, drop = FALSE]
    }
  }
  if(!is.null(classes))
  list(dataTable = dataTable[, -match(c("primary", "class"), colnames(dataTable))],
       classes = classes)
  else
    dataTable[, -match("primary", colnames(dataTable))]
}

.checkVariablesAndSame <- function(trainingMatrix, testingMatrix)
{
  if(ncol(trainingMatrix) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else if(ncol(trainingMatrix) != ncol(testingMatrix))
    stop("Training dataset and testing dataset contain differing numbers of features.")  
}

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
      if("trainParams" %in% names(.methodFormals(selectParams@featureSelection))) # Needs training and prediction functions for resubstitution error rate calculation.
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
        paramList <- append(paramList, c(selParams, datasetName = "N/A", selectionName = "N/A"))
        do.call(selector, paramList)
      }, selectParams@featureSelection, selectParams@otherParams, SIMPLIFY = FALSE)
      
      if(class(featuresLists[[1]]) == "SelectResult")
      {
        featuresCounts <- table(unlist(lapply(featuresLists, function(featureSet) featureSet@chosenFeatures[[1]])))
        selectedFeatures <- as.integer(names(featuresCounts))[featuresCounts >= selectParams@minPresence]
      } else { # The prediction function used for resubstitution returned a variety of lists.
        selectedFeatures <- lapply(1:length(featuresLists[[1]]), function(variety)
        {
          varietyFeatures <- lapply(featuresLists, function(selectList) selectList[[variety]]@chosenFeatures[[1]])
          featuresCounts <- table(unlist(varietyFeatures))
          as.integer(names(featuresCounts))[featuresCounts >= selectParams@minPresence]
        })
        names(selectedFeatures) <- names(featuresLists[[1]])
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

.doTrain <- function(measurements, classes, training, testing, trainParams,
                            predictParams, verbose)
  # Re-use inside feature selection.
{
  initialClass <- class(measurements)
  if(class(expression) != "list") # Will be a DataFrame.
    measurements <- list(data = measurements)

  trained <- mapply(function(measurementsVariety, variety)
  {
    measurementsTrain <- measurementsVariety[training, ]
    measurementsTest <- measurementsVariety[testing, ]
    
    if(variety != "data") # Single expression set is in a list with name 'data'.
    {
      multiplierParams <- sapply(strsplit(variety, ",")[[1]], strsplit, split = '=')
      individiualParams <- lapply(multiplierParams, '[', 2)
      names(individiualParams) <- sapply(multiplierParams, '[', 1)
      individiualParams <- lapply(individiualParams, function(param) tryCatch(as.numeric(param), warning = function(warn){param}))
      trainFormals <- tryCatch(names(.methodFormals(trainParams@classifier@generic)), warning = function(warn) {names(formals(trainParams@classifier))})
      predictFormals <- tryCatch(names(formals(predictParams@predictor)), warning = function(warn) {tryCatch(names(.methodFormals(predictParams@predictor)), character(0))}) 
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
    paramList <- list(measurementsTrain, classes[training])
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
            paramList <- list(model, measurementsTrain) # Model and test set same as training set.
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
              calcExternalPerformance(classes[training], predictions, resubstituteParams@performanceType)
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
      paramList <- append(paramList, list(measurementsTest))
      if(length(predictParams@otherParams) > 0)
      {
        useOthers <- setdiff(names(predictParams@otherParams), names(trainParams@otherParams))
        paramList <- c(paramList, trainParams@otherParams, predictParams@otherParams[useOthers])
      }
      paramList <- c(paramList, verbose = verbose)
            
      if(is.null(tuneCombinations))
      {
        predictions <- do.call(trainParams@classifier, paramList)
        if(verbose >= 2)
          message("Training and prediction completed.")    
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
            performanceValue <- calcExternalPerformance(classes[training], predictedFactor[[predictIndex]],
                                                        resubstituteParams@performanceType)
            list(predictedClasses[[predictIndex]], performanceValue)
          })
        })
        
        chosenPredictions <- lapply(1:length(performances[[1]]), function(predictVariety)
        {
          performanceValues <- sapply(performances, function(tuneLevel) tuneLevel[[predictVariety]][[2]]) # Value is in second position.
          if(resubstituteParams@better == "lower")
            chosenTune <- which.min(performanceValues)[1]
          else
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

.pickFeatures <- function(measurements, classes, datasetName, trainParams, predictParams,
                          resubstituteParams, orderedFeatures, selectionName, verbose)
{
  performances <- sapply(resubstituteParams@nFeatures, function(topFeatures)
  {
    measurementsSubset <- measurements[, orderedFeatures[1:topFeatures]]
    trained <- .doTrain(measurementsSubset, classes, 1:nrow(measurementsSubset), 1:nrow(measurementsSubset),
                        trainParams, predictParams, verbose)
    
    if(trainParams@doesTests == FALSE)
    {
        predictions <- .doTest(trained, measurementsSubset, 1:nrow(measurementsSubset),
                               predictParams, verbose)
    } else { # Same function does training and testing.
      if(class(trained) != "list" || sum(grepl('=', names(trained))) == 0)
        predictions <- predictParams@getClasses(trained)
      else
        predictions <- lapply(trained, function(trainedModel) predictParams@getClasses(trainedModel))
    }
    if(class(predictions) == "list") # Mutiple varieties of predictions.
    {
      if(class(predictions[[1]]) == "data.frame")
        labels <- lapply(predictions, function(set) set[, sapply(set, class) == "factor"])
      else
        labels <- predictions
    } else { # A single variety of prediction.
      if(class(predictions) == "data.frame")
        labels <- predictions[, sapply(predictions, class) == "factor"]
      else
        labels <- predictions
    }

    if(class(labels) == "list")
    {
      performanceValues <- lapply(labels, function(labelSet) calcExternalPerformance(classes, labelSet, resubstituteParams@performanceType))
      
    } else {
      performanceValues <- calcExternalPerformance(classes, labels, resubstituteParams@performanceType)
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
  
  rankedFeatures <- lapply(1:length(pickedFeatures), function(variety) orderedFeatures)
  pickedFeatures <- lapply(pickedFeatures, function(pickedSet) orderedFeatures[pickedSet])
  
  if(!is.null(mcols(measurements))) # Table describing source table and variable name is present.
  {
    varInfo <- mcols(measurements)
    rankedFeatures <- lapply(rankedFeatures, function(features) varInfo[features, ])
    pickedFeatures <- lapply(pickedFeatures, function(features) varInfo[features, ])
  }
  
  selectResults <- lapply(1:length(rankedFeatures), function(variety)
  {
    SelectResult(datasetName, selectionName,
                 list(rankedFeatures[[variety]]), list(pickedFeatures[[variety]]))
  })
  names(selectResults) <- names(pickedFeatures)
  
  if(length(selectResults) == 1) selectResults <- selectResults[[1]] else selectResults
}

.validationText <- function(result)
{
  switch(result@validation[[1]],
  resampleFold = paste(result@validation[[2]], "Permutations,", result@validation[[3]], "Folds"),
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

.methodFormals <- function(f, signature = "ExpressionSet") {
  tryCatch({
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
    genFormals},
    error = function(error) {
      formals(f)
    })
}

.densityCrossover <- function(aDensity, anotherDensity)
{
  if(!all(aDensity[['x']] == anotherDensity[['x']]))
    stop("x positions are not the same for the two density variables.")
  allDifferences <- aDensity[['y']] - anotherDensity[['y']]
  crosses <- which(diff(sign(allDifferences)) != 0)
  if(aDensity[['y']][crosses[1]] == 0 && aDensity[['y']][crosses[length(crosses)]] == 0)
    crosses <- crosses[-c(1, length(crosses))] # Remove crossings at ends of densities.
  aDensity[['x']][crosses]
}
