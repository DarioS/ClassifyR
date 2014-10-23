setOldClass("pamrtrained")

.doSelection <- function(expression, training, selectionParams, trainParams,
                         predictParams, verbose)
{
  initialClass <- class(expression)
  if(class(expression) != "list")
    expression <- list(data = expression)  

  selected <- lapply(expression, function(expressionVariety)
  {
    if(is.function(selectionParams@featureSelection))
    {
      paramList <- list(expressionVariety[, training], trainParams = trainParams,
                        predictParams = predictParams, verbose = verbose)
      paramList <- append(paramList, selectionParams@otherParams)
      do.call(selectionParams@featureSelection, paramList)    
    } else { # It is a list of functions.
      featuresLists <- mapply(function(selector, selParams)
      {
        paramList <- list(expressionVariety[, training], trainParams = trainParams,
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
  })

  if(initialClass != "list") selected <- selected[[1]]
  if(class(selected[[1]]) == "list") selected <- unlist(selected, recursive = FALSE)
  selected
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
    
    if(variety != "data")
    {
      multiplierParams <- sapply(strsplit(variety, ",")[[1]], strsplit, split = '=')
      individiualParams <- lapply(multiplierParams, '[', 2)
      names(individiualParams) <- sapply(multiplierParams, '[', 1)
      individiualParams <- lapply(individiualParams, function(param) tryCatch(as.numeric(param), warning = function(warn){param}))
      if(all(names(individiualParams) %in% names(trainParams@otherParams)))
        trainParams@otherParams[names(individiualParams)] <- individiualParams
      else
        predictParams@otherParams[names(individiualParams)] <- individiualParams
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
      trained
    } else { # Some classifiers do training and testing with a single function.
      paramList <- append(paramList, list(expressionTest))
      if(length(predictParams@otherParams) > 0)
        paramList <- c(paramList, predictParams@otherParams)
      paramList <- c(paramList, verbose = verbose)
      prediction <- do.call(trainParams@classifier, paramList)
      if(verbose >= 2)
        message("Training and classification completed.")    
      prediction
    }
  }, expression, names(expression), SIMPLIFY = FALSE)
  
  if(initialClass != "list") trained <- trained[[1]]
  if(class(trained[[1]]) == "list") trained <- unlist(trained, recursive = FALSE)
  trained
}

.doTest <- function(trained, expression, testing, predictParams, verbose)
# Re-use inside feature selection.
{
  initialClass <- class(trained)
  if(class(trained) != "list")
  {
    trained <- list(model = trained)
    expression <- list(data = expression)
  }
  
  predicted <- mapply(function(model, data)
  {
    if(!grepl("{}", paste(capture.output(predictParams@predictor), collapse = ''), fixed = TRUE))
    {
      testExpression <- exprs(data)[, testing]
      if(predictParams@transposeExpression == TRUE)
        testExpression <- t(testExpression)
    
      paramList <- list(model, testExpression)
      if(length(predictParams@otherParams) > 0)
        paramList <- c(paramList, predictParams@otherParams)
      paramList <- c(paramList, verbose = verbose)
       prediction <- do.call(predictParams@predictor, paramList)
    } else
    {
      prediction <- model
    }
    predictions <- predictParams@getClasses(prediction)
    
    if(verbose >= 2)
      message("Prediction completed.")    
    predictions
  }, trained, expression, SIMPLIFY = FALSE)

  if(initialClass != "list") predicted <- predicted[[1]]
  if(class(predicted[[1]]) == "list") predicted <- unlist(predicted, recursive = FALSE)
  predicted
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
