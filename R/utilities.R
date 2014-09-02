.doSelection <- function(expression, training, selectionParams, trainParams,
                         predictParams, verbose)
{
  if(is.function(selectionParams@featureSelection))
  {
    paramList <- list(expression[, training], trainParams = trainParams,
                      predictParams = predictParams, verbose = verbose)
    paramList <- append(paramList, selectionParams@otherParams)
    do.call(selectionParams@featureSelection, paramList)    
  } else { # It is a list of functions.
    featuresLists <- mapply(function(selector, selParams)
    {
      paramList <- list(expression[, training], trainParams = trainParams,
                        predictParams = predictParams, verbose = verbose)
      paramList <- append(paramList, selParams)
      do.call(selector, paramList)    
    }, selectionParams@featureSelection, selectionParams@otherParams, SIMPLIFY = FALSE)
    if(is.numeric(featuresLists[[1]]))
    {
      featuresCounts <- table(unlist(featuresLists))
      as.integer(names(featuresCounts))[featuresCounts >= selectionParams@minPresence]
    } else {
      lapply(1:length(featuresLists[[1]]), function(variety)
      {
        varietyFeatures <- table(unlist(lapply(featuresLists, "[[", variety)))
        as.integer(names(varietyFeatures))[varietyFeatures >= selectionParams@minPresence]
      })
    }
  }
}

.doTransform <- function(expression, transformParams, verbose)
{
  if(!grepl("{}", paste(capture.output(transformParams@transform), collapse = ''), fixed = TRUE))
  {
    paramList <- list(expression)
    do.call(transformParams@transform,
            c(paramList, transformParams@otherParams, verbose = verbose))
  } else {
    expression # Return expression, unchanged.
  }  
}

.doTrainAndTest <- function(expression, training, testing, trainParams,
                            predictParams, verbose)
# Re-use inside feature selection.
{
  initialClass <- class(expression)
  if(!is.list(expression))
    expression <- list(data = expression)

  varPredictions <- mapply(function(expressionVariety, variety)
  {
    classes <- pData(expressionVariety)[, "class"]
    expressionTrain <- exprs(expressionVariety[, training])
    expressionTest <- exprs(expressionVariety[, testing])
    if(trainParams@transposeExpression == TRUE) # For classifiers that consider columns as samples and genes as rows.
    {
      expressionTrain <- t(expressionTrain)
      expressionTest <- t(expressionTest)
    }
  
    paramList <- list(expressionTrain, classes[training])
    if(trainParams@doesTests == FALSE) # Training and prediction are separate.
    {
      if(length(trainParams@otherParams) > 0)
        paramList <- c(paramList, trainParams@otherParams)    
      paramList <- c(paramList, verbose = verbose)
      trained <- do.call(trainParams@classifier, paramList)
      if(verbose >= 2)
        message("Training completed.")  
      paramList <- list(trained, expressionTest)
      if(length(predictParams@otherParams) > 0)
         paramList <- c(paramList, predictParams@otherParams)
      paramList <- c(paramList, verbose = verbose)
      prediction <- do.call(predictParams@predictor, paramList)
      if(verbose >= 2)
        message("Prediction completed.")
    } else { # Some classifiers do training and testing with a single function.
      paramList <- append(paramList, list(expressionTest))
      if(length(predictParams@otherParams) > 0)
        paramList <- c(paramList, predictParams@otherParams)
      paramList <- c(paramList, verbose = verbose)
      prediction <- do.call(trainParams@classifier, paramList)
      if(verbose >= 2)
        message("Classification completed.")    
    }
    predictions <- predictParams@getClasses(prediction)
    if(initialClass == "list") predictions[[variety]] else predictions
  }, expression, names(expression), SIMPLIFY = FALSE)
  if(initialClass == "ExpressionSet") varPredictions[[1]] else varPredictions
}

.pickRows <- function(errorRates)
{
  if(is.numeric(errorRates))
  {
    nFeatures <- as.numeric(names(errorRates))
    pickedRows <- 1:(nFeatures[which.min(errorRates)[1]])
  } else {
    nFeatures <- as.numeric(colnames(errorRates))
    pickedRows <- apply(errorRates, 1, function(rates) 1:(nFeatures[which.min(rates)[1]]))
    if(is.matrix(pickedRows))
      pickedRows <- as.list(as.data.frame(pickedRows))
  }
  pickedRows
}