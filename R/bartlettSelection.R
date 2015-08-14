setGeneric("bartlettSelection", function(expression, ...)
{standardGeneric("bartlettSelection")})

setMethod("bartlettSelection", "matrix", 
          function(expression, classes, ...)
          {
            colnames(expression) <- NULL # Might be duplicates because of sampling with replacement. 
            features <- rownames(expression)
            groupsTable <- data.frame(class = classes)
            exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
            if(length(features) > 0) featureNames(exprSet) <- features  
            bartlettSelection(ExpressionSet(expression, AnnotatedDataFrame(groupsTable)), ...)
          })

setMethod("bartlettSelection", "ExpressionSet", 
          function(expression, datasetName, trainParams, predictParams, resubstituteParams, selectionName = "Bartlett Test", verbose = 3)
          {
            if(verbose == 3)
              message("Calculating Bartlett statistic.")
            exprMatrix <- exprs(expression)
            classes <- pData(expression)[, "class"]
            pValues <- apply(exprMatrix, 1, function(geneRow) bartlett.test(geneRow, classes)[["p.value"]])
            orderedFeatures <- order(pValues)
            
            .pickRows(expression, datasetName, trainParams, predictParams, resubstituteParams, orderedFeatures, selectionName, verbose)
          })