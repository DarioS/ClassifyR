setGeneric("leveneSelection", function(expression, ...)
           {standardGeneric("leveneSelection")})

setMethod("leveneSelection", "matrix", 
          function(expression, classes, ...)
{
  colnames(expression) <- NULL # Might be duplicates because of sampling with replacement. 
  features <- rownames(expression)
  groupsTable <- data.frame(class = classes)
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features  
  leveneSelection(ExpressionSet(expression, AnnotatedDataFrame(groupsTable)), ...)
})

setMethod("leveneSelection", "ExpressionSet", 
          function(expression, nFeatures, trainParams, predictParams, verbose = 3)
{
  if(!requireNamespace("car", quietly = TRUE))
    stop("The package 'car' could not be found. Please install it.")            
  if(verbose == 3)
    message("Calculating Levene statistic.")
  exprMatrix <- exprs(expression)
  classes <- pData(expression)[, "class"]
  pValues <- apply(exprMatrix, 1, function(geneRow) car::leveneTest(geneRow, classes)[["Pr(>F)"]][1])
  orderedFeatures <- order(pValues)
  if(verbose == 3)
    message("Selecting number of features to use.")
  errorRates <- sapply(nFeatures, function(topFeatures)
  {
    expressionSubset <- expression[orderedFeatures[1:topFeatures], ]
    trained <- .doTrain(expressionSubset, 1:ncol(expressionSubset), 1:ncol(expressionSubset),
                        trainParams, predictParams, verbose)
    if(trainParams@doesTests == FALSE)
      predictions <- .doTest(trained, expressionSubset, 1:ncol(expressionSubset),
                             predictParams, verbose)
    else
      predictions <- predictParams@getClasses(trained)
    
    if(is.list(predictions))
      lapply(predictions, function(predictions) sum(predictions != classes) / length(classes))
    else
      sum(predictions != classes) / length(classes)
  })
  if(class(errorRates) == "numeric") names(errorRates) <- nFeatures else colnames(errorRates) <- nFeatures
  
  picked <- .pickRows(errorRates)
  if(verbose == 3)
    message("Features selected.")
  
  if(class(picked) == "list")
  {
    rankedFeatures <- lapply(1:length(picked), function(variety) orderedFeatures)
    pickedFeatures <- lapply(picked, function(pickedSet) orderedFeatures[pickedSet])
  } else {
    rankedFeatures <- orderedFeatures
    pickedFeatures <- orderedFeatures[picked]
  }
  list(rankedFeatures, pickedFeatures)
})
