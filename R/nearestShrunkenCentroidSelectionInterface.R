setGeneric("nearestShrunkenCentroidSelectionInterface", function(expression, ...)
{standardGeneric("nearestShrunkenCentroidSelectionInterface")})

setMethod("nearestShrunkenCentroidSelectionInterface", "matrix", function(expression, classes, ...)
{
  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")

  colnames(expression) <- NULL # Might be duplicates because of sampling with replacement.
  features <- rownames(expression)
  groupsTable <- data.frame(class = classes)
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features
  nearestShrunkenCentroidSelectionInterface(exprSet, ...)
})

setMethod("nearestShrunkenCentroidSelectionInterface", "ExpressionSet",
          function(expression, datasetName, trained, ..., selectionName = "Shrunken Centroids",
                   verbose = 3)
{ 
  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")
  
  classes <- pData(expression)[, "class"]
  minError <- min(trained[["errors"]])
  threshold <- trained[["threshold"]][max(which(trained[["errors"]] == minError))]
  
  params <- list(...)
  params <- params[!names(params) %in% c("trainParams", "predictParams")]
  params <- c(list(trained), list(list(x = exprs(expression), y = classes, geneid = 1:nrow(expression))), threshold, params)
  
  chosen <- as.numeric(do.call(pamr::pamr.listgenes, params)[, 1])
  
  if(verbose == 3)
    message("Nearest shrunken centroid feature selection completed.")
  
  SelectResult(datasetName, selectionName, list(), list(chosen))
})