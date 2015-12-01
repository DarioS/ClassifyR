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
                   absolute = TRUE, verbose = 3)
{
  location <- match.arg(location)
  expressionTrain <- exprs(expression)[, training]
  if(location == "mean")
    geneTrainingLocations <- rowMeans(expressionTrain)
  else # median.
    geneTrainingLocations <- apply(expressionTrain, 1, median)
  transformed <- apply(exprs(expression), 2, '-', geneTrainingLocations)
  if(absolute == TRUE)
    transformed <- abs(transformed)
  exprs(expression) <- transformed
  if(verbose == 3)
    message("Subtraction from ", location,
            {if(absolute == TRUE) " and absolute transformation"}, " completed.")
  
  expression
})