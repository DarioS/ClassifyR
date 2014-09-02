setGeneric("KullbackLeiblerSelection", function(expression, ...)
           {standardGeneric("KullbackLeiblerSelection")})

setMethod("KullbackLeiblerSelection", "matrix", function(expression, classes, ...)
{ 
  colnames(expression) <- NULL # Might be duplicates because of sampling with replacement.  
  features <- rownames(expression)
  rownames(expression) <- NULL
  groupsTable <- data.frame(class = classes)
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features
  KullbackLeiblerSelection(exprSet, ...)
})

setMethod("KullbackLeiblerSelection", "ExpressionSet", 
          function(expression, nFeatures, trainParams,
                                     predictParams, type = c("symmetric", "dynamic"),
                                     ..., verbose = 3)
{
  type <- match.arg(type)
  if(verbose == 3)
    message("Selecting features by Kullback Leibler divergence")
  classes <- pData(expression)[, "class"]
  oneClassExpression <- exprs(expression[, classes == levels(classes)[1]])
  otherClassExpression <- exprs(expression[, classes == levels(classes)[2]])
  oneClassDistribution <- getLocationsAndScales(oneClassExpression, ...)
  otherClassDistribution <- getLocationsAndScales(otherClassExpression, ...)
  locationDifference <- oneClassDistribution[[1]] - otherClassDistribution[[1]]
  if(type == "symmetric")
    divergence <- 1/2 * (locationDifference^2/((oneClassDistribution[[2]])^2) +
                         locationDifference^2/((otherClassDistribution[[2]])^2) +
                         ((oneClassDistribution[[2]])^2) / ((otherClassDistribution[[2]])^2) +
                         ((otherClassDistribution[[2]])^2) / ((oneClassDistribution[[2]])^2))
  else
    divergence <- abs(locationDifference) + abs(oneClassDistribution[[2]]
                - otherClassDistribution[[2]])
  
  orderedFeatures <- order(divergence, decreasing = TRUE)
  if(verbose == 3)
    message("Selecting number of features to use.")

  errorRates <- sapply(nFeatures, function(topFeatures)
  {
    expressionSubset <- expression[orderedFeatures[1:topFeatures], ]
    trainPredictions <- .doTrainAndTest(expressionSubset, 1:ncol(expressionSubset), 1:ncol(expressionSubset),
                                     trainParams, predictParams, verbose = verbose)
    if(is.list(trainPredictions))
      lapply(trainPredictions, function(predictions) sum(predictions != classes) / length(classes))
    else
      sum(trainPredictions != classes) / length(classes)
  })
  if(class(errorRates) == "numeric") names(errorRates) <- nFeatures else colnames(errorRates) <- nFeatures
  
  picked <- .pickRows(errorRates)
  if(verbose == 3)
    message("Features selected.")

  if(class(picked) == "list")
    lapply(picked, function(pickedSet) orderedFeatures[pickedSet])
  else
    orderedFeatures[picked]
})