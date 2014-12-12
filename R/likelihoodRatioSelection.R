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
          function(expression, trainParams, predictParams, resubstituteParams,
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
  
  .pickRows(expression, trainParams, predictParams, resubstituteParams, orderedFeatures, verbose)
})
