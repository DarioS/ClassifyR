setGeneric("subtractFromLocation", function(expression, ...)
           {standardGeneric("subtractFromLocation")})

setMethod("subtractFromLocation", "matrix", 
          function(expression, ...)
{ 
  colnames(expression) <- NULL
  rownames(expression) <- NULL
  subtractFromLocation(ExpressionSet(expression), ...)
})

setMethod("subtractFromLocation", "ExpressionSet", 
          function(expression, training, location = c("mean", "median"),
                   verbose = 3)
{
  location <- match.arg(location)
  expressionTrain <- exprs(expression)[, training]
  if(location == "mean")
    geneTrainingLocations <- rowMeans(expressionTrain)
  else # median.
    geneTrainingLocations <- apply(expressionTrain, 1, median)
  transformed <- apply(expression, 2, '-', geneTrainingLocations)
  exprs(expression) <- transformed
  if(verbose == 3)
    message("Subtraction from ", location, " completed.")
  
  expression
})