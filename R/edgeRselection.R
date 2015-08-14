setGeneric("edgeRselection", function(expression, ...)
           {standardGeneric("edgeRselection")})

setMethod("edgeRselection", "matrix", 
          function(expression, classes, ...)
{
  colnames(expression) <- NULL # Might be duplicates because of sampling with replacement.            
  features <- rownames(expression)
  groupsTable <- data.frame(class = classes)
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features
  edgeRselection(exprSet, ...)
})

setMethod("edgeRselection", "ExpressionSet", 
          function(expression, datasetName, normFactorsOptions = NULL,
                   dispOptions = NULL, fitOptions = NULL, trainParams,
                   predictParams, resubstituteParams, selectionName = "edgeR LRT",
                   verbose = 3)
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

  .pickRows(expression, datasetName, trainParams, predictParams, resubstituteParams, orderedFeatures,
            selectionName, verbose)
})
