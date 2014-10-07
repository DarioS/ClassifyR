setGeneric("limmaSelection", function(expression, ...)
           {standardGeneric("limmaSelection")})

setMethod("limmaSelection", "matrix", 
          function(expression, classes, ...)
{
  colnames(expression) <- NULL # Might be duplicates because of sampling with replacement.   
  features <- rownames(expression)
  groupsTable <- data.frame(class = classes)   
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features
  limmaSelection(ExpressionSet(expression, AnnotatedDataFrame(groupsTable)), ...)
})

setMethod("limmaSelection", "ExpressionSet", 
          function(expression, nFeatures, trainParams, predictParams, ..., verbose = 3)
{
  if(!requireNamespace("limma", quietly = TRUE))
    stop("The package 'limma' could not be found. Please install it.")            
  if(verbose == 3)
    message("Doing feature selection.")
  exprMatrix <- exprs(expression)
  classes <- pData(expression)[, "class"]
  allFeatures <- featureNames(expression)
  
  fitParams <- list(exprMatrix, model.matrix(~ classes))
  if(!missing(...))
    fitParams <- append(fitParams, list(...))
  prognosisModel <- do.call(limma::lmFit, fitParams)
  prognosisModel <- limma::eBayes(prognosisModel)
  orderedFeatures <- match(rownames(limma::topTable(prognosisModel, 2, number = Inf, sort.by = "p")), allFeatures)
  
  if(verbose == 3)
    message("Selecting number of features to use.")
  errorRates <- sapply(nFeatures, function(topFeatures)
  {
    expressionSubset <- expression[orderedFeatures[1:topFeatures], ]
    trained <- .doTrain(expressionSubset, 1:ncol(expressionSubset), 1:ncol(expressionSubset),
                        trainParams, predictParams, verbose)
    if(trainParams@doesTests == FALSE)
      predictions <- .doTest(trained, expressionSubset, 1:ncol(expressionSubset),
                             predictParams, verbose)
    else
      predictions <- trained
    
    if(is.list(predictions))
      lapply(predictions, function(predictions) sum(predictions != classes) / length(classes))
    else
      sum(predictions != classes) / length(classes)
  })
  if(class(errorRates) == "numeric") names(errorRates) <- nFeatures else colnames(errorRates) <- nFeatures
  
  picked <- .pickRows(errorRates)
  if(verbose == 3)
    message("Features selected.")
  
  if(class(picked) == "list")
    lapply(picked, function(pickedSet) orderedFeatures[pickedSet])
  else
    orderedFeatures[picked]
})