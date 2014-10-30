setGeneric("likelihoodRatioSelection", function(expression, ...)
           {standardGeneric("likelihoodRatioSelection")})

setMethod("likelihoodRatioSelection", "matrix", function(expression, classes, ...)
{ 
  colnames(expression) <- NULL # Might be duplicates because of sampling with replacement.  
  features <- rownames(expression)
  groupsTable <- data.frame(class = classes)
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features
  likelihoodRatioSelection(exprSet, ...)
})

setMethod("likelihoodRatioSelection", "ExpressionSet", 
          function(expression, nFeatures, trainParams, predictParams,
                   alternative = c(location = "different", scale = "different"),
                   ..., verbose = 3)
{
  if(verbose == 3)
    message("Selecting features by likelihood ratio ranking.")
  
  classes <- pData(expression)[, "class"]
  exprMatrix <- exprs(expression)
  oneClass <- classes == levels(classes)[1]
  otherClass <- classes == levels(classes)[2]
  oneClassExpression <- exprMatrix[, oneClass]
  otherClassExpression <- exprMatrix[, otherClass]
  oneClassDistribution <- getLocationsAndScales(oneClassExpression, ...)
  otherClassDistribution <- getLocationsAndScales(otherClassExpression, ...)
  allDistribution <- getLocationsAndScales(exprMatrix, ...)

  logLikelihoodRatios <- -2 * (unlist(mapply(function(geneRow, scale, location)
  sum(dnorm(geneRow, scale, location, log = TRUE)),
  as.data.frame(t(exprMatrix)), allDistribution[[1]], allDistribution[[2]])) -
  unlist(mapply(function(geneRow, scale, location)
  sum(dnorm(geneRow, scale, location, log = TRUE)),
  as.data.frame(t(oneClassExpression)),
  switch(alternative[["location"]], same = allDistribution[[1]], different = oneClassDistribution[[1]]),
  switch(alternative[["scale"]], same = allDistribution[[2]], different = oneClassDistribution[[2]]))) +
  unlist(mapply(function(geneRow, scale, location)
  sum(dnorm(geneRow, scale, location, log = TRUE)),
  as.data.frame(t(otherClassExpression)),
  switch(alternative[["location"]], same = allDistribution[[1]], different = otherClassDistribution[[1]]),
  switch(alternative[["scale"]], same = allDistribution[[2]], different = otherClassDistribution[[2]]))))

  orderedFeatures <- order(logLikelihoodRatios)
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
