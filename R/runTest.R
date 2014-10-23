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
                               if(length(transformParams@intermediate) != 0)
                                 transformParams@otherParams <- c(transformParams@otherParams,
                                                                  as.list(environment())[transformParams@intermediate])

                               transformParams@otherParams <- c(transformParams@otherParams, list(training = training))
                               expression <- .doTransform(expression, transformParams, verbose)
                               newSize <- if(class(expression) == "list") length(expression) else 1
                             },
                 SelectionParams = {
                               if(length(selectionParams@intermediate) != 0)
                                 selectionParams@otherParams <- c(selectionParams@otherParams,
                                                                  as.list(environment())[selectionParams@intermediate])
                               
                               selectedFeatures <- .doSelection(expression, training, selectionParams,
                                                                trainParams, predictParams, verbose)

                               if(selectionParams@subsetExpressionData == TRUE)
                               {
                                 if(is.numeric(selectedFeatures) || is.character(selectedFeatures))
                                 {
                                   if(class(expression) != "list")
                                     expression <- expression[selectedFeatures, ]
                                   else
                                     lapply(expression, function(variety) variety[selectedFeatures, ])
                                 } else {
                                   if(class(expression) != "list")
                                     expression <- lapply(selectedFeatures, function(features) expression[features, ])
                                   else
                                     expression <- lapply(expression, function(variety)
                                                          lapply(selectedFeatures, function(features) variety[features, ]))
                                 }
                               } else {
                                 if(is.list(selectedFeatures))
                                 {
                                   if(class(expression) != "list")
                                     expression <- lapply(selectedFeatures, function(features) expression)
                                   else
                                     expression <- lapply(expression, function(variety)
                                                          lapply(selectedFeatures, function(features) variety))
                                 }
                               }
                               
                               if(class(expression) == "list" && class(expression[[1]]) == "list")
                                 expression <- unlist(expression, recursive = FALSE)
                               newSize <- length(expression)
                               lastSize <- newSize
                             }, 
                 TrainParams = {
                              if(length(trainParams@intermediate) != 0)
                                trainParams@otherParams <- c(trainParams@otherParams,
                                                             as.list(environment())[trainParams@intermediate])

                              trained <- .doTrain(expression, training, testing, trainParams, predictParams, verbose)
                              newSize <- if(class(trained) == "list") length(trained) else 1
                              if(newSize / lastSize != 1) expression <- unlist(lapply(expression, function(variety)
                                                                                      lapply(1:(newSize / lastSize), function(x) variety)),
                                                                               recursive = FALSE)
                              lastSize <- newSize
                              },
                 PredictParams = {
                                 if(length(predictParams@intermediate) != 0)
                                   predictParams@otherParams <- c(predictParams@otherParams,
                                                      as.list(environment())[predictParams@intermediate])
                                  predictedClasses <- .doTest(trained, expression, testing, predictParams, verbose)
                                 }
           )
    
  }
  
  list(selectedFeatures, testing, predictedClasses)
})
