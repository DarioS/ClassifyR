setGeneric("leveneSelection", function(expression, ...)
           {standardGeneric("leveneSelection")})

setMethod("leveneSelection", "matrix", 
          function(expression, classes, ...)
{
  colnames(expression) <- NULL # Might be duplicates because of sampling with replacement. 
  features <- rownames(expression)
  rownames(expression) <- NULL # In case of duplicate gene symbols rownames.
  groupsTable <- data.frame(class = classes)
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features  
  leveneSelection(ExpressionSet(expression, AnnotatedDataFrame(groupsTable)), ...)
})

setMethod("leveneSelection", "ExpressionSet", 
          function(expression, nFeatures, trainParams, predictParams, verbose = 3)
{
  if(verbose == 3)
    message("Calculating Levene statistic.")
  exprMatrix <- exprs(expression)
  classes <- pData(expression)[, "class"]
  pValues <- apply(exprMatrix, 1, function(geneRow) leveneTest(geneRow, classes)[["Pr(>F)"]][1])
  orderedFeatures <- order(pValues)
  if(verbose == 3)
    message("Selecting number of features to use.")
  errorRates <- sapply(nFeatures, function(topFeatures)
  {
    expressionSubset <- expression[orderedFeatures[1:topFeatures], ]
    sum(.doTrainAndTest(expressionSubset, 1:ncol(expressionSubset), 1:ncol(expressionSubset),
                        trainParams, predictParams, verbose = verbose) != classes) / length(classes)
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