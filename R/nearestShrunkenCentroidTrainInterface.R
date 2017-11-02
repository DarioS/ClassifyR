setGeneric("nearestShrunkenCentroidTrainInterface", function(expression, ...)
{standardGeneric("nearestShrunkenCentroidTrainInterface")})

setMethod("nearestShrunkenCentroidTrainInterface", "matrix", function(expression, classes, ...)
{
  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")
  
  features <- rownames(expression)
  groupsTable <- data.frame(class = classes, row.names = colnames(expression))
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features
  nearestShrunkenCentroidTrainInterface(exprSet, ...)
})

setMethod("nearestShrunkenCentroidTrainInterface", "ExpressionSet", function(expression, ..., verbose = 3)
{ 
  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")
  
  classes <- pData(expression)[, "class"]
  trainedModel <- pamr::pamr.train(list(x = exprs(expression), y = classes,
                                        geneid = 1:nrow(expression)), ...)
  
  if(verbose == 3)
    message("Nearest shrunken centroid training completed.")
  
  trainedModel
})