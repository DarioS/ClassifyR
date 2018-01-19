setGeneric("edgeRselection", function(counts, ...)
           {standardGeneric("edgeRselection")})

setMethod("edgeRselection", "matrix", # Matrix of integer counts.
          function(counts, classes, ...)
{
  if(is.null(names(classes))) names(classes) <- colnames(counts)
  edgeRselection(MultiAssayExperiment(list(matrix = counts), DataFrame(class = classes)),
                 targets = "matrix", ...)
})

# One or more omics datasets, possibly with clinical data.
setMethod("edgeRselection", "MultiAssayExperiment",
          function(counts, datasetName, targets = NULL,
                   normFactorsOptions = NULL, dispOptions = NULL, fitOptions = NULL,
                   trainParams, predictParams, resubstituteParams,
                   selectionName = "edgeR LRT", verbose = 3)
{
  if(!requireNamespace("edgeR", quietly = TRUE))
    stop("The package 'edgeR' could not be found. Please install it.")
  if(is.null(targets))
    stop("'targets' must be specified but was not.")
  if(length(setdiff(targets, names(counts))))
    stop("Some values of 'targets' are not names of 'counts' but all must be.")            
            
  if(verbose == 3)
    message("Doing feature selection.")

  tablesAndClasses <- .MAEtoWideTable(counts, targets, "integer")
  countsTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]

  # DGEList stores features as rows and samples as columns.          
  countsList <- edgeR::DGEList(t(as.matrix(countsTable)), group = classes)
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
  orderedFeatures <- match(rownames(result[["table"]]), colnames(countsTable))

  .pickFeatures(countsTable, classes, datasetName,
                trainParams, predictParams, resubstituteParams,
                orderedFeatures, selectionName, verbose)
})
