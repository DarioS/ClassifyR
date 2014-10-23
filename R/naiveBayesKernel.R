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
                   minDifference = 0, verbose = 3)
{
  weighted <- match.arg(weighted)
            
  classes <- pData(expression)[, "class"]
  classesSizes <- sapply(levels(classes), function(class) sum(classes == class))
  expression <- exprs(expression)      
  
  if(verbose == 3)
    message("Fitting densities.")
  classesSplines <- lapply(levels(classes), function(class)
  {
    classDensities <- apply(expression[, classes == class], 1, function(geneRow) 
                      {
                        geneDensity <- density(geneRow)
                        splinefun(geneDensity[['x']], geneDensity[['y']], "natural")
                      })
  })
  
  if(verbose == 3)
    message("Predicting classes using fitted densities.")
  classesDensities <- lapply(classesSplines, function(classes)
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
    sapply(posterior, function(sampleCol)
    {
      sampleCol <- sampleCol[abs(sampleCol) > difference]
      sum(sampleCol > 0) > length(sampleCol) / 2
    })
  })
  
  weightedOther <- lapply(minDifference, function(difference)
  {
    sapply(posterior, function(sampleCol)
    {
      sampleCol <- sampleCol[abs(sampleCol) > difference]
      sum(sampleCol) > 0
    })
  })
  
  unweightedList <- lapply(lapply(unweightedOther, function(other)
                    {  
                      predictions <- rep(levels(classes)[1], ncol(test))
                      predictions[other] <- levels(classes)[2]
                      predictions
                    }), factor, levels = levels(classes))
  weightedList <- lapply(lapply(weightedOther, function(other)
                  {  
                    predictions <- rep(levels(classes)[1], ncol(test))
                    predictions[other] <- levels(classes)[2]
                    predictions  
                  }), factor, levels = levels(classes))

  if(length(unweightedList) == 1)
  {
    unweightedList <- unlist(unweightedList)
    weightedList <- unlist(weightedList)
  } else {
    names(unweightedList) <- paste("weighted=unweighted,minDifference=", minDifference, sep = '')
    names(weightedList) <- paste("weighted=weighted,minDifference=", minDifference, sep = '')
  }
                              
  switch(weighted, unweighted = unweightedList,
                   weighted = weightedList,
                   both = if(class(unweightedList) == "list")
                            unlist(list(unweightedList, weightedList), recursive = FALSE)
                          else
                            list(Unweighted = unweightedList, Weighted = weightedList)
                       
        )
})
