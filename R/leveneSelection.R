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
          function(expression, trainParams, predictParams, resubstituteParams, verbose = 3)
{
  if(!requireNamespace("car", quietly = TRUE))
    stop("The package 'car' could not be found. Please install it.")            
  if(verbose == 3)
    message("Calculating Levene statistic.")
  exprMatrix <- exprs(expression)
  classes <- pData(expression)[, "class"]
  pValues <- apply(exprMatrix, 1, function(geneRow) car::leveneTest(geneRow, classes)[["Pr(>F)"]][1])
  orderedFeatures <- order(pValues)
 
  .pickRows(expression, trainParams, predictParams, resubstituteParams, orderedFeatures, verbose)
})
