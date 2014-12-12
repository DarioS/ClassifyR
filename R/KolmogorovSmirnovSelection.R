setGeneric("KolmogorovSmirnovSelection", function(expression, ...)
           {standardGeneric("KolmogorovSmirnovSelection")})

setMethod("KolmogorovSmirnovSelection", "matrix", function(expression, classes, ...)
{ 
  colnames(expression) <- NULL # Might be duplicates because of sampling with replacement.  
  features <- rownames(expression)
  groupsTable <- data.frame(class = classes)
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features
  KolmogorovSmirnovSelection(exprSet, ...)
})

setMethod("KolmogorovSmirnovSelection", "ExpressionSet", 
          function(expression, trainParams, predictParams, resubstituteParams, ..., verbose = 3)
{
  if(verbose == 3)
    message("Selecting features by Kolmogorov Smirnov distance")
  classes <- pData(expression)[, "class"]
  oneClass <- classes == levels(classes)[1]
  otherClass <- classes == levels(classes)[2]
  KSdistance <- apply(exprs(expression), 1, function(geneRow)
                      ks.test(geneRow[oneClass], geneRow[otherClass], ...)[["statistic"]])

  orderedFeatures <- order(KSdistance, decreasing = TRUE)
  .pickRows(expression, trainParams, predictParams, resubstituteParams, orderedFeatures, verbose)
})