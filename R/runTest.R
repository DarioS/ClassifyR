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
          function(expression, training, testing,
                   params = list(SelectionParams(), TrainParams(), PredictParams()),
                   verbose = 1)
{
  stagesParamClasses <- sapply(params, class)
  if(match("TrainParams", stagesParamClasses) > match("PredictParams", stagesParamClasses))
    stop("\"testing\" must not be before \"training\" in 'params'.")
  
  transformParams <- params[[match("TransformParams", stagesParamClasses)]]
  selectionParams <- params[[match("SelectionParams", stagesParamClasses)]]
  trainParams <- params[[match("TrainParams", stagesParamClasses)]]
  predictParams <- params[[match("PredictParams", stagesParamClasses)]]
  
  lastSize <- 1
  for(stageIndex in 1:length(params))
  {
    switch(stagesParamClasses[[stageIndex]],
                 TransformParams = {
                               transformParams@otherParams <- c(transformParams@otherParams, list(training = training))
                               expression <- .doTransform(expression, transformParams, verbose)
                               newSize <- if(class(expression) == "list") length(expression) else 1
                             },
                 SelectionParams = {
                               selectedFeatures <- .doSelection(expression, training, selectionParams,
                                                                trainParams, predictParams, verbose)
                               if(is.numeric(selectedFeatures) || is.character(selectedFeatures))
                                 expression <- expression[selectedFeatures, ]
                               else
                                 expression <- lapply(selectedFeatures, function(features) expression[features, ])
                               newSize <- if(class(selectedFeatures) == "list") length(selectedFeatures) else 1
                               if(newSize / lastSize != 1) expression <- unlist(lapply(expression, function(variety)
                                                                                       lapply(1:(newSize / lastSize), function(x) variety)),
                                                                                recursive = FALSE)                               
                               lastSize <- newSize
                             }, 
                 TrainParams = {
                              trained <- .doTrain(expression, training, testing, trainParams, predictParams, verbose)
                              newSize <- if(class(trained) == "list") length(trained) else 1
                              if(newSize / lastSize != 1) expression <- unlist(lapply(expression, function(variety)
                                                                                      lapply(1:(newSize / lastSize), function(x) variety)),
                                                                               recursive = FALSE)
                              lastSize <- newSize
                              },
                 PredictParams = predictedClasses <- .doTest(trained, expression, testing, predictParams, verbose)
           )
    
  }
  
  list(selectedFeatures, testing, predictedClasses)
})