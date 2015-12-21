setGeneric("previousSelection", function(expression, ...)
{standardGeneric("previousSelection")})

setMethod("previousSelection", "matrix", 
          function(expression, classes, ...)
          {
            colnames(expression) <- NULL # Might be duplicates because of sampling with replacement. 
            features <- rownames(expression)
            groupsTable <- data.frame(class = classes)
            exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
            if(length(features) > 0) featureNames(exprSet) <- features  
            previousSelection(ExpressionSet(expression, AnnotatedDataFrame(groupsTable)), ...)
          })

setMethod("previousSelection", "ExpressionSet", 
          function(expression, datasetName, classifyResult, minimumOverlapPercent = 80,
                   selectionName = "Previous Selection", .iteration, verbose = 3)
          {
            if(verbose == 3)
              message("Choosing previous features.")
            
            if(length(.iteration) == 1)
              previousIDs <- featureNames(classifyResult)[features(classifyResult)[[.iteration]]]
            else # Resample index and fold index.
              previousIDs <- featureNames(classifyResult)[features(classifyResult)[[.iteration[[1]]]][[.iteration[[2]]]]]
            
            indicesInCurrent <- match(previousIDs, rownames(expression))
            commonFeatures <- sum(!is.na(indicesInCurrent)) / length(indicesInCurrent) * 100
            if(commonFeatures < minimumOverlapPercent)
              signalCondition(simpleError("Number of features in common between previous and current dataset is lower than 'minimumOverlapPercent'."))
            
            SelectResult(datasetName, selectionName, list(), list(na.omit(indicesInCurrent))) # Ranking isn't transferred across.
          })