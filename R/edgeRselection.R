setGeneric("edgeRselection", function(expression, ...)
           {standardGeneric("edgeRselection")})

setMethod("edgeRselection", "matrix", 
          function(expression, classes, ...)
{
  colnames(expression) <- NULL # Might be duplicates because of sampling with replacement.            
  features <- rownames(expression)
  rownames(expression) <- NULL # In case of duplicate gene symbols rownames.
  groupsTable <- data.frame(class = classes)
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features
  edgeRselection(exprSet, ...)
})

setMethod("edgeRselection", "ExpressionSet", 
          function(expression, nFeatures, normFactorsOptions = NULL,
                           dispOptions = NULL, fitOptions = NULL, trainParams,
                           predictParams, verbose = 3)
{
  if(!requireNamespace("edgeR", quietly = TRUE))
    stop("The package 'edgeR' could not be found. Please install it.")
  if(verbose == 3)
    message("Doing feature selection.")
  
  allFeatures <- featureNames(expression)
  exprMatrix <- exprs(expression)
  colnames(exprMatrix) <- NULL # Might be duplicates because of sampling with replacement.  
  classes <- pData(expression)[, "class"]
  expressionList <- edgeR::DGEList(exprMatrix, group = classes)
  paramList <- list(expressionList)
  if(!is.null(normFactorsOptions))
    paramList <- append(paramList, normFactorsOptions)
  if(verbose == 3)
    message("Calculating scaling factors.")
  expressionList <- do.call(edgeR::calcNormFactors, paramList)
  paramList <- list(expressionList, model.matrix(~ classes))
  if(!is.null(dispOptions))
    paramList <- append(paramList, dispOptions)
  if(verbose == 3)
    message("Estimating dispersion.")
  expressionList <- do.call(edgeR::estimateDisp, paramList)
  paramList <- list(expressionList, model.matrix(~ classes))
  if(!is.null(fitOptions))
    paramList <- append(paramList, fitOptions)
  if(verbose == 3)
    message("Fitting linear model.")
  fit <- do.call(edgeR::glmFit, paramList)
  result <- edgeR::topTags(edgeR::glmLRT(fit, coef = 2), n = Inf, adjust.method = "none")
  orderedFeatures <- match(rownames(result[["table"]]), allFeatures)

  if(verbose == 3)
    message("Selecting number of features to use.")
  errorRates <- sapply(nFeatures, function(topFeatures)
  {
    expressionSubset <- expression[orderedFeatures[1:topFeatures], ]
    sum(.doTrainAndTest(expressionSubset, 1:ncol(expressionSubset), 1:ncol(expressionSubset),
                      trainParams, predictParams, verbose = verbose) != classes) / length(classes)
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