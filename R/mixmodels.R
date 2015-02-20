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
          function(models, test, weighted = c("both", "unweighted", "weighted"),
                   minDifference = 0, returnType = c("label", "score", "both"), verbose = 3)
{
  weighted <- match.arg(weighted)
  returnType <- match.arg(returnType)
  test <- exprs(test)
  classLevels <- attr(models, "classes")
  
  if(verbose == 3)
    message("Predicting using normal mixtures.")
  
  otherClassScore <- mapply(function(geneExpr, oneClassModel, otherClassModel)
                     {
                       classScores <- lapply(list(oneClassModel, otherClassModel), function(model)
                                      {
                                        Reduce('+', lapply(1:model@bestResult@nbCluster, function(index)
                                        {
                                          model@bestResult@parameters@proportions[index] * dnorm(geneExpr, model@bestResult@parameters@mean[index], sqrt(as.numeric(model@bestResult@parameters@variance[[index]])))
                                        }))
                                      })
                       models[[2]][[2]] * classScores[[2]] - models[[1]][[2]] * classScores[[1]]
                       # Second element of second list in 'models' is unimportant information added by mixmod.
                }, data.frame(t(test)), models[[1]][[1]], models[[2]][[1]]) 
  
  otherClassScore <- as.data.frame(otherClassScore)

  unweightedOther <- lapply(minDifference, function(difference)
  {
    apply(otherClassScore, 1, function(sampleRow)
    {
      sampleRow <- sampleRow[abs(sampleRow) > difference]
      list(sum(sampleRow > 0) > length(sampleRow) / 2, sum(sampleRow > 0) / length(sampleRow))
    })
  })
  unweightedIsOther <- lapply(unweightedOther, function(difference) unlist(lapply(difference, "[[", 1)))
  unweightedOtherScores <- lapply(unweightedOther, function(difference) unlist(lapply(difference, "[[", 2)))

  weightedOther <- lapply(minDifference, function(difference)
  {
    weightedOtherClass <- apply(otherClassScore, 1, function(sampleRow)
    {
      sampleRow <- sampleRow[abs(sampleRow) > difference]
      list(sum(sampleRow) > 0, sum(sampleRow))
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
