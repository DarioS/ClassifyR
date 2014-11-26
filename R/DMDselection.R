setGeneric("DMDselection", function(expression, ...)
           {standardGeneric("DMDselection")})

setMethod("DMDselection", "matrix", function(expression, classes, ...)
{ 
  colnames(expression) <- NULL # Might be duplicates because of sampling with replacement.  
  features <- rownames(expression)
  groupsTable <- data.frame(class = classes)
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features
  DMDselection(exprSet, ...)
})

setMethod("DMDselection", "ExpressionSet", 
          function(expression, nFeatures, trainParams,
                                     predictParams, ..., verbose = 3)
{
  if(verbose == 3)
    message("Selecting features by DMD.")
  classes <- pData(expression)[, "class"]
  oneClassExpression <- exprs(expression[, classes == levels(classes)[1]])
  otherClassExpression <- exprs(expression[, classes == levels(classes)[2]])
  oneClassDistribution <- getLocationsAndScales(oneClassExpression, ...)
  otherClassDistribution <- getLocationsAndScales(otherClassExpression, ...)
  locationDifference <- oneClassDistribution[[1]] - otherClassDistribution[[1]]
  divergence <- abs(locationDifference) + abs(oneClassDistribution[[2]] - otherClassDistribution[[2]])
  
  orderedFeatures <- order(divergence, decreasing = TRUE)
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
      predictions <- trained
    
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