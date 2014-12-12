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
          function(expression, trainParams, predictParams, resubstituteParams, ..., verbose = 3)
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
  
  .pickRows(expression, trainParams, predictParams, resubstituteParams, orderedFeatures, verbose)
})