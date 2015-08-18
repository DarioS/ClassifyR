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
  
  other <- lapply(minDifference, function(difference)
  {
    lapply(posterior, function(sampleCol)
    {
      sampleCol <- sampleCol[abs(sampleCol) > difference]
      if(length(sampleCol) == 0) # No features have a large enough density difference.
      {                          # Simply vote for the larger class.
        largerClass <- names(classesSizes)[which.max(classesSizes)]
        if(largerClass == levels(classes)[1])
        {
          logicalSymbol <- FALSE
          difference <- -1
        } else {
          logicalSymbol <- TRUE
          difference <- 1
        }
        list(unweighted=list(logicalSymbol, difference),
             weighted=list(logicalSymbol, difference))
      } else { 
        list(unweighted=list(sum(sampleCol > 0) > length(sampleCol) / 2, sum(sampleCol > 0) / length(sampleCol)),
             weighted=list(sum(sampleCol) > 0, sum(sampleCol)))
      }
    })
  })
  
  unweightedIsOther <- lapply(other, function(difference) unlist(lapply(difference, function(sample) sample[["unweighted"]][[1]])))
  unweightedOtherScores <- lapply(other, function(difference) unlist(lapply(difference, function(sample) sample[["unweighted"]][[2]])))
  weightedIsOther <- lapply(other, function(difference) unlist(lapply(difference, function(sample) sample[["weighted"]][[1]])))
  weightedOtherScores <- lapply(other, function(difference) unlist(lapply(difference, function(sample) sample[["weighted"]][[2]])))
         
  predictionsList <- mapply(function(weightVarietyLabels, weightVarietyScores)
                     {  
                       mapply(function(other, scores)
                       {
                       predictions <- rep(levels(classes)[1], ncol(test))
                       predictions[other] <- levels(classes)[2]
                       predictions <- factor(predictions, levels = levels(classes))
                       predictions
                       switch(returnType, label = predictions, score = scores,
                                          both = data.frame(label = predictions, score = scores))
                       }, weightVarietyLabels, weightVarietyScores, SIMPLIFY = FALSE)
                     }, list(unweightedIsOther, weightedIsOther),
                        list(unweightedOtherScores, weightedOtherScores), SIMPLIFY = FALSE)

  if(length(predictionsList[[1]]) == 1) # No minDifference range.
  {
    if(class(predictionsList[[1]][[1]]) != "data.frame")
    {
      predictionsList <- lapply(predictionsList, unlist)
    } else {
      predictionsList <- lapply(predictionsList, "[[", 1)
    }
  } else {
    names(predictionsList[[1]]) <- paste("weighted=unweighted,minDifference=", minDifference, sep = '')
    names(predictionsList[[2]]) <- paste("weighted=weighted,minDifference=", minDifference, sep = '')
  }
  
  switch(weighted, unweighted = predictionsList[[1]],
         weighted = predictionsList[[2]],
         both = if(class(predictionsList[[1]]) == "list")
           unlist(list(predictionsList[[1]], predictionsList[[2]]), recursive = FALSE)
         else
           list(`weighted=unweighted` = predictionsList[[1]], `weighted=weighted` = predictionsList[[2]])
         )
})