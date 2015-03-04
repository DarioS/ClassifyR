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
          function(expression, datasetName, classificationName, training, testing,
                   params = list(SelectionParams(), TrainParams(), PredictParams()),
                   verbose = 1, .usedInternally = FALSE)
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
                               
                               topFeatures <- .doSelection(expression, training, selectionParams,
                                                                trainParams, predictParams, verbose)

                               if(sum(grepl('=', names(topFeatures[[1]]))) > 0)
                               {
                                 multiSelection <- TRUE
                               } else {
                                 multiSelection <- FALSE
                               }
                               
                               rankedFeatures <- topFeatures[[1]] # Extract for adding to result list.
                               selectedFeatures <- topFeatures[[2]]

                               if(selectionParams@subsetExpressionData == TRUE)
                               {
                                 if(multiSelection == FALSE)
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
                                 if(is.list(selectedFeatures)) # Multiple selection varieties. Replicate the expression data.
                                 {
                                   if(class(expression) != "list")
                                     expression <- lapply(selectedFeatures, function(features) expression)
                                   else
                                     expression <- lapply(expression, function(variety)
                                                          lapply(selectedFeatures, function(features) variety))
                                 }
                               }

                               if(class(expression) == "list" && class(expression[[1]]) == "list")
                               {
                                 oldNames <- sapply(expression, names)
                                 newNames <- unlist(lapply(expression, names))
                                 expression <- unlist(expression, recursive = FALSE)
                                 names(expression) <- paste(rep(oldNames, each = length(expression[[1]])), newNames, sep = ',')
                               }
                               newSize <- length(expression)
                               lastSize <- newSize
                             }, 
                 TrainParams = {
                              if(length(trainParams@intermediate) != 0)
                                trainParams@otherParams <- c(trainParams@otherParams,
                                                             as.list(environment())[trainParams@intermediate])
                              
                              trained <- .doTrain(expression, training, testing, trainParams, predictParams, verbose)

                              newSize <- if(class(trained) == "list") length(trained) else 1
                              if(newSize / lastSize != 1)
                              {
                                expression <- unlist(lapply(if(class(expression) == "list") expression else list(expression), function(variety)
                                                                                      lapply(1:(newSize / lastSize), function(replicate) variety)),
                                                                               recursive = FALSE)
                                names(expression) <- names(trained)
                              }
                              
                              lastSize <- newSize
                              if(class(trained) == "list")
                                tuneDetails <- lapply(trained, attr, "tune")
                              else
                                tuneDetails <- attr(trained, "tune")
                                if(is.null(tuneDetails)) tuneDetails <- list(tuneDetails)
                              },
                 PredictParams = {
                                 if(length(predictParams@intermediate) != 0)
                                   predictParams@otherParams <- c(predictParams@otherParams,
                                                      as.list(environment())[predictParams@intermediate])
                                  predictedClasses <- .doTest(trained, expression, testing, predictParams, verbose)
                                 }
           )
    
  }
  if(class(testing) == "logical") testing <- which(testing)

  if(.usedInternally == TRUE)
  {        
    list(rankedFeatures, selectedFeatures, testing, predictedClasses, tuneDetails)
  } else { # runTest is being used directly, rather than from runTests. Create a ClassifyResult object.
    if(class(predictedClasses) != "list")
    {
      return(ClassifyResult(datasetName, classificationName, sampleNames(expression), featureNames(expression),
                            list(rankedFeatures), list(selectedFeatures), list(data.frame(sample = testing, label = predictedClasses)),
                            pData(expression)[, "class"], list("independent"), tuneDetails)
             )
    } else {
      return(mapply(function(varietyPredictions, varietyTunes)
      {
        if(is.null(varietyTunes)) varietyTunes <- list(varietyTunes)
        ClassifyResult(datasetName, classificationName, sampleNames(expression), featureNames(expression),
                       list(rankedFeatures), list(selectedFeatures), list(data.frame(sample = testing, label = varietyPredictions)),
                       pData(expression)[, "class"], list("independent"), varietyTunes)
      }, predictedClasses, tuneDetails, SIMPLIFY = FALSE))
    }
  }
})
