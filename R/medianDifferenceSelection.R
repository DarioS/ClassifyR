setGeneric("medianDifferenceSelection", function(expression, ...)
{standardGeneric("medianDifferenceSelection")})

setMethod("medianDifferenceSelection", "matrix", 
          function(expression, classes, ...)
          {
            colnames(expression) <- NULL # Might be duplicates because of sampling with replacement. 
            features <- rownames(expression)
            groupsTable <- data.frame(class = classes)
            exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
            if(length(features) > 0) featureNames(exprSet) <- features  
            medianDifferenceSelection(ExpressionSet(expression, AnnotatedDataFrame(groupsTable)), ...)
          })

setMethod("medianDifferenceSelection", "ExpressionSet", 
          function(expression, datasetName, trainParams, predictParams, resubstituteParams, selectionName = "Difference of Group Medians", verbose = 3)
          {
            if(verbose == 3)
              message("Calculating differences of group medians.")
            exprMatrix <- exprs(expression)
            classes <- pData(expression)[, "class"]
            classExprMatrix <- exprMatrix[, classes == levels(classes)[1]]
            otherClassExprMatrix <- exprMatrix[, classes == levels(classes)[2]]
            medianAbsoluteDifferences <- abs(apply(classExprMatrix, 1, median) - apply(classExprMatrix, 2, median))
            orderedFeatures <- order(medianAbsoluteDifferences, decreasing = TRUE)
            
            .pickRows(expression, datasetName, trainParams, predictParams, resubstituteParams, orderedFeatures, selectionName, verbose)
          })