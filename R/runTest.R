setGeneric("runTest", function(expression, ...)
           {standardGeneric("runTest")})

setMethod("runTest", c("matrix"),
  function(expression, classes, ...)
{
    groupsTable <- data.frame(class = classes)
    rownames(groupsTable) <- colnames(expression)
    exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
    runTest(exprSet, ...)
})
setMethod("runTest", c("ExpressionSet"),
          function(expression, doFirst = c("transform", "selection"), training, testing,
                   transformParams = TransformParams(), selectionParams = SelectionParams(),
                   trainParams = TrainParams(), predictParams = PredictParams(), verbose = 1)
{    
  doFirst <- match.arg(doFirst)
  if(!grepl("{}", paste(capture.output(transformParams@transform), collapse = ''), fixed = TRUE))
    transformParams@otherParams <- c(transformParams@otherParams, list(training = training))
  if(doFirst == "transform") # Transformation, followed by feature selection.
  {
    transformed <- .doTransform(expression, transformParams, verbose)
    selectedFeatures <- .doSelection(transformed, training, selectionParams,
                                     trainParams, predictParams, verbose)
    if(is.numeric(selectedFeatures) || is.character(selectedFeatures))
      expression <- transformed[selectedFeatures, ]
    else
      expression <- lapply(selectedFeatures, function(features) transformed[features, ])
  } else { # Feature selection, then transformation.
    selectedFeatures <- .doSelection(expression, training, selectionParams,
                                     trainParams, predictParams, verbose)

    expression <- .doTransform(expression, transformParams, verbose)
    if(is.numeric(selectedFeatures) || is.character(selectedFeatures))  
      expression <- expression[selectedFeatures, ]
    else
      expression <- lapply(selectedFeatures, function(features) expression[features, ])
  }

  predictedClasses <- .doTrainAndTest(expression, training, testing,
                                      trainParams, predictParams, verbose)  
  list(selectedFeatures, testing, predictedClasses)
})