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
          function(expression, trainParams, predictParams, resubstituteParams, ..., verbose = 3)
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
    fitParams <- append(fitParams, ...)
  prognosisModel <- do.call(limma::lmFit, fitParams)
  prognosisModel <- limma::eBayes(prognosisModel)
  orderedFeatures <- match(rownames(limma::topTable(prognosisModel, 2, number = Inf, sort.by = "p")), allFeatures)
  
  .pickRows(expression, trainParams, predictParams, resubstituteParams, orderedFeatures, verbose)
})
