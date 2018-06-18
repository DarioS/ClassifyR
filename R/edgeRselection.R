setGeneric("edgeRselection", function(counts, ...)
           {standardGeneric("edgeRselection")})

setMethod("edgeRselection", "matrix", # Matrix of integer counts.
          function(counts, classes, ...)
{
  edgeRselection(DataFrame(t(counts), check.names = FALSE), classes, ...)
})

# DataFrame of counts, likely created by runTests or runTest.
setMethod("edgeRselection", "DataFrame",
          function(counts, classes, datasetName,
                   normFactorsOptions = NULL, dispOptions = NULL, fitOptions = NULL,
                   trainParams, predictParams, resubstituteParams,
                   selectionName = "edgeR LRT", verbose = 3)
{
  if(verbose == 3)
    message("Doing edgeR LRT feature selection.")
  
  # DGEList stores features as rows and samples as columns.          
  countsList <- edgeR::DGEList(t(as.matrix(counts)), group = classes)
  paramList <- list(countsList)
  if(!is.null(normFactorsOptions))
    paramList <- append(paramList, normFactorsOptions)
  if(verbose == 3)
    message("Calculating scaling factors.")
  countsList <- do.call(edgeR::calcNormFactors, paramList)
  paramList <- list(countsList, model.matrix(~ classes))
  if(!is.null(dispOptions))
    paramList <- append(paramList, dispOptions)
  if(verbose == 3)
    message("Estimating dispersion.")
  countsList <- do.call(edgeR::estimateDisp, paramList)
  paramList <- list(countsList, model.matrix(~ classes))
  if(!is.null(fitOptions))
    paramList <- append(paramList, fitOptions)
  if(verbose == 3)
    message("Fitting linear model.")
  fit <- do.call(edgeR::glmFit, paramList)
  result <- edgeR::topTags(edgeR::glmLRT(fit, coef = 2), n = Inf, adjust.method = "none")
  orderedFeatures <- match(rownames(result[["table"]]), colnames(counts))
  
  .pickFeatures(counts, classes, NULL, datasetName,
                trainParams, predictParams, resubstituteParams,
                orderedFeatures, selectionName, verbose)    
})

# One or more omics data sets, possibly with clinical data.
setMethod("edgeRselection", "MultiAssayExperiment",
          function(counts, targets = NULL, ...)
{
  if(!requireNamespace("edgeR", quietly = TRUE))
    stop("The package 'edgeR' could not be found. Please install it.")
  if(is.null(targets))
    stop("'targets' must be specified but was not.")
  if(length(setdiff(targets, names(counts))))
    stop("Some values of 'targets' are not names of 'counts' but all must be.")            

  tablesAndClasses <- .MAEtoWideTable(counts, targets, "integer")
  countsTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  edgeRselection(countsTable, classes, ...)
})