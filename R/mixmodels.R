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
                   minDifference = 0, verbose = 3)
{
  weighted <- match.arg(weighted)            
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
                }, data.frame(t(test)), models[[1]][[1]], models[[2]][[1]])
  
  otherClassScore <- as.data.frame(otherClassScore)

  unweightedOther <- lapply(minDifference, function(difference)
  {
    apply(otherClassScore, 1, function(sampleRow)
    {
      sampleRow <- sampleRow[abs(sampleRow) > difference]
      sum(sampleRow > 0) > length(sampleRow) / 2
    })
  })
  
  weightedOther <- lapply(minDifference, function(difference)
  {
    weightedOtherClass <- apply(otherClassScore, 1, function(sampleRow)
    {
      sampleRow <- sampleRow[abs(sampleRow) > difference]
      sum(sampleRow) > 0
    })
  })
  
  unweightedList <- lapply(lapply(unweightedOther, function(other)
  {  
    predictions <- rep(classLevels[1], ncol(test))
    predictions[other] <- classLevels[2]
    predictions
  }), factor, levels = classLevels)
  weightedList <- lapply(lapply(weightedOther, function(other)
  {  
    predictions <- rep(classLevels[1], ncol(test))
    predictions[other] <- classLevels[2]
    predictions  
  }), factor, levels = classLevels)
  
  if(length(unweightedList) == 1)
  {
    unweightedList <- unlist(unweightedList)
    weightedList <- unlist(weightedList)
  } else {
    names(unweightedList) <- paste("Unweighted,minDifference", minDifference, sep = '=')
    names(weightedList) <- paste("Weighted,minDifference", minDifference, sep = '=')
  }
  
  switch(weighted, unweighted = unweightedList,
         weighted = weightedList,
         both = if(class(unweightedList) == "list")
           unlist(list(unweightedList, weightedList), recursive = FALSE)
         else
           list(Unweighted = unweightedList, Weighted = weightedList)
         
  )
})