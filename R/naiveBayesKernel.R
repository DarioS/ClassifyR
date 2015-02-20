setGeneric("naiveBayesKernel", function(expression, ...)
           {standardGeneric("naiveBayesKernel")})

setMethod("naiveBayesKernel", "matrix", 
          function(expression, classes, ...)
{ 
  colnames(expression) <- NULL # Might be duplicates because of sampling with replacement.            
  features <- rownames(expression)
  groupsTable <- data.frame(class = classes)  
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features
  naiveBayesKernel(exprSet, ...)
})

setMethod("naiveBayesKernel", "ExpressionSet", 
          function(expression, test, weighted = c("both", "unweighted", "weighted"),
                   minDifference = 0, returnType = c("label", "score", "both"), verbose = 3)
{
  weighted <- match.arg(weighted)
  returnType <- match.arg(returnType)
            
  classes <- pData(expression)[, "class"]
  classesSizes <- sapply(levels(classes), function(class) sum(classes == class))
  expression <- exprs(expression)      
  
  if(verbose == 3)
    message("Fitting densities.")
  classesSplines <- lapply(levels(classes), function(class) # Two lists for each class of
                                                            # densities for each feature.
  {
    classDensities <- apply(expression[, classes == class], 1, function(geneRow) 
                      {
                        geneDensity <- density(geneRow)
                        splinefun(geneDensity[['x']], geneDensity[['y']], "natural")
                      })
  })
  
  if(verbose == 3)
    message("Predicting classes using fitted densities.")
  classesDensities <- lapply(classesSplines, function(classes) # Density values at test features' expression.
  {
    densities <- as.matrix(mapply(function(geneSpline, geneSamples)
    {
      geneSpline(geneSamples)
    }, classes, as.data.frame(t(test))))
    if(ncol(test) > 1) densities <- t(densities) else densities
  })
  
  posterior <- classesSizes[2] * classesDensities[[2]] - classesSizes[1] * classesDensities[[1]]
  posterior <- as.data.frame(posterior)
  
  unweightedOther <- lapply(minDifference, function(difference)
  {
    as.data.frame(sapply(posterior, function(sampleCol)
    {
      sampleCol <- sampleCol[abs(sampleCol) > difference]
      list(sum(sampleCol > 0) > length(sampleCol) / 2, sum(sampleCol > 0) / length(sampleCol))
    }))
  })
  unweightedIsOther <- lapply(unweightedOther, function(difference) unlist(lapply(difference, "[[", 1)))
  unweightedOtherScores <- lapply(unweightedOther, function(difference) unlist(lapply(difference, "[[", 2)))

  weightedOther <- lapply(minDifference, function(difference)
  {
    lapply(posterior, function(sampleCol)
    {
      sampleCol <- sampleCol[abs(sampleCol) > difference]
      list(sum(sampleCol) > 0, sum(sampleCol))
    })
  })
  weightedIsOther <- lapply(weightedOther, function(difference) unlist(lapply(difference, "[[", 1)))
  weightedOtherScores <- lapply(weightedOther, function(difference) unlist(lapply(difference, "[[", 2)))
                                
  unweightedList <- mapply(function(other, scores)
                    {  
                      predictions <- rep(levels(classes)[1], ncol(test))
                      predictions[other] <- levels(classes)[2]
                      predictions <- factor(predictions, levels = levels(classes))
                      predictions
                      switch(returnType, label = predictions, score = scores,
                                         both = data.frame(label = predictions, score = scores))
                    }, unweightedIsOther, unweightedOtherScores, SIMPLIFY = FALSE)
  weightedList <- mapply(function(other, scores)
                  {  
                    predictions <- rep(levels(classes)[1], ncol(test))
                    predictions[other] <- levels(classes)[2]
                    predictions <- factor(predictions, levels = levels(classes))
                    predictions
                    switch(returnType, label = predictions, score = scores,
                           both = data.frame(label = predictions, score = scores))
                  }, weightedIsOther, weightedOtherScores, SIMPLIFY = FALSE)

  if(length(unweightedList) == 1)
  {
    if(class(unweightedList[[1]]) != "data.frame")
    {
      unweightedList <- unlist(unweightedList)
      weightedList <- unlist(weightedList)
    } else {
      unweightedList <- unweightedList[[1]]
      weightedList <- weightedList[[1]]
    }
  } else {
    names(unweightedList) <- paste("weighted=unweighted,minDifference=", minDifference, sep = '')
    names(weightedList) <- paste("weighted=weighted,minDifference=", minDifference, sep = '')
  }
  
  switch(weighted, unweighted = unweightedList,
         weighted = weightedList,
         both = if(class(unweightedList) == "list")
           unlist(list(unweightedList, weightedList), recursive = FALSE)
         else
           list(`weighted=unweighted` = unweightedList, `weighted=weighted` = weightedList)
         
  )
})