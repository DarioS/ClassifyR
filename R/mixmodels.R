setGeneric("mixModelsTrain", function(expression, ...)
           {standardGeneric("mixModelsTrain")})

setMethod("mixModelsTrain", "matrix", function(expression, classes, ...)
{ 
  colnames(expression) <- NULL # Might be duplicates because of sampling with replacement.
  features <- rownames(expression)
  groupsTable <- data.frame(class = classes)
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features
  mixModelsTrain(exprSet, ...)
})

setMethod("mixModelsTrain", "ExpressionSet",
          function(expression, ..., verbose = 3)
{
  if(verbose == 3)
    message("Fitting mixtures of normals for genes.")
  if(!requireNamespace("Rmixmod", quietly = TRUE))
    stop("The package 'Rmixmod' could not be found. Please install it.")
  classes <- pData(expression)[, "class"]
  expression <- exprs(expression)
  oneClass <- classes == levels(classes)[1]
  otherClass <- classes == levels(classes)[2]
  oneClassExpression <- expression[, oneClass]
  otherClassExpression <- expression[, otherClass]
  isExtra <- !missing(...)
  
  models <- lapply(list(oneClassExpression, otherClassExpression), function(classExpression)
            {
                classExpression <- data.frame(t(classExpression))
                list(apply(classExpression, 2, function(geneRow)
                {
                   mixmodParams <- list(geneRow)
                   if(isExtra) mixmodParams <- append(mixmodParams, list(...))
                   do.call(Rmixmod::mixmodCluster, mixmodParams)
                }), nrow(classExpression))        
            })
  
  if(verbose == 3)
    message("Done fitting normal mixtures.")
  attr(models, "classes") <- levels(classes)
  models
})

setGeneric("mixModelsTest", function(models, test, ...)
           {standardGeneric("mixModelsTest")})

setMethod("mixModelsTest", c("list", "matrix"), function(models, test, ...)
{ 
  colnames(test) <- NULL # Might be duplicates because of sampling with replacement.  
  exprSet <- ExpressionSet(test)
  mixModelsTest(models, exprSet, ...)
})

setMethod("mixModelsTest", c("list", "ExpressionSet"),
          function(models, test, densityXvalues = 1024,
                   weighted = c("both", "unweighted", "weighted"),
                   weight = c("both", "height difference", "crossover distance"),
                   minDifference = 0, tolerance = 0.01,
                   returnType = c("label", "score", "both"), verbose = 3)
{
  weighted <- match.arg(weighted)
  returnType <- match.arg(returnType)
  test <- exprs(test)
  classLevels <- attr(models, "classes")
  
  if(verbose == 3)
    message("Predicting using normal mixtures.")
  
  otherClassScore <- mapply(function(featureValues, oneClassModel, otherClassModel)
                     {
    browser()
                       featureValues <- c(oneClassModel@data, otherClassModel@data)
                       xValues <- seq(min(featureValues), max(featureValues), length.out = densityXvalues)
                       # Convert to a density, so the crossover points can be calculated.
                       classScores <- lapply(list(oneClassModel, otherClassModel), function(model)
                                      {
                                        Reduce('+', lapply(1:model@bestResult@nbCluster, function(index)
                                        {
                                          model@bestResult@parameters@proportions[index] * dnorm(featureValues, model@bestResult@parameters@mean[index], sqrt(as.numeric(model@bestResult@parameters@variance[[index]])))
                                        }))
                                      })
                       models[[2]][[2]] * classScores[[2]] - models[[1]][[2]] * classScores[[1]]
                       # Second element of second list in 'models' is unimportant information added by mixmod.
                }, data.frame(t(test)), models[[1]][[1]], models[[2]][[1]]) 
  
  otherClassScore <- as.data.frame(otherClassScore)

  other <- lapply(minDifference, function(difference)
  {
    apply(otherClassScore, 1, function(sampleRow)
    {
      sampleRow <- sampleRow[abs(sampleRow) > difference]
      if(length(sampleRow) == 0) # No features have a large enough density difference.
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
        list(unweighted=list(sum(sampleRow > 0) > length(sampleRow) / 2, sum(sampleRow > 0) / length(sampleRow)),
             weighted=list(sum(sampleRow) > 0, sum(sampleRow)))
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
