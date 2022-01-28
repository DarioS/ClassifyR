#' Reproducibly Run Various Kinds of Cross-Validation
#' 
#' Enables doing classification schemes such as ordinary 10-fold, 100
#' permutations 5-fold, and leave one out cross-validation.  Processing in
#' parallel is possible by leveraging the package \code{\link{BiocParallel}}.
#' 
#' 
#' @aliases runTests runTests,matrix-method runTests,DataFrame-method
#' runTests,MultiAssayExperiment-method
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix}, the rows are features, and the columns are samples.  The
#' sample identifiers must be present as column names of the \code{matrix} or
#' the row names of the \code{DataFrame}. 
#' @param classes A vector of class labels of class \code{\link{factor}} of the
#' same length as the number of samples in \code{measurements} if it is a
#' \code{\link{matrix}} (i.e. number of columns) or a \code{\link{DataFrame}}
#' (i.e. number of rows) or a character vector of length 1 containing the
#' column name in \code{measurements} if it is a \code{\link{DataFrame}} or the
#' column name in \code{colData(measurements)} if \code{measurements} is a
#' \code{\link{MultiAssayExperiment}}. If a column name, that column will be
#' removed before training.
#' @param crossValParams An object of class \code{\link{CrossValParams}},
#' specifying the kind of cross-validation to be done.
#' @param modellingParams An object of class \code{\link{ModellingParams}},
#' specifying the class rebalancing, transformation (if any), feature selection
#' (if any), training and prediction to be done on the data set.
#' @param characteristics A \code{\link{DataFrame}} describing the
#' characteristics of the classification used. First column must be named
#' \code{"charateristic"} and second column must be named \code{"value"}.
#' Useful for automated plot annotation by plotting functions within this
#' package.  Transformation, selection and prediction functions provided by
#' this package will cause the characteristics to be automatically determined
#' and this can be left blank.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"clinical"} is also a valid value
#' and specifies that the clinical data table will be used.
#' @param outcomeColumn The column name of the clinical data table containing
#' the outcome to be predicted. Will automatically be removed from the clinical
#' data table during model training. Must be specified if data is a
#' \code{MultiAssayExperiment}.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method.
#' @param verbose Default: 1. A number between 0 and 3 for the amount of
#' progress messages to give.  A higher number will produce more messages as
#' more lower-level functions print messages.
#' @return An object of class \code{\link{ClassifyResult}}.
#' @author Dario Strbenac
#' @examples
#' 
#'   #if(require(sparsediscrim))
#'   #{
#'     data(asthma)
#'     
#'     CVparams <- CrossValParams(permutations = 5)
#'     tuneList <- list(nFeatures = seq(5, 25, 5), performanceType = "Balanced Error")
#'     selectParams <- SelectParams(differentMeansRanking, tuneParams = tuneList)
#'     modellingParams <- ModellingParams(selectParams = selectParams)
#'     runTests(measurements, classes, CVparams, modellingParams,
#'              DataFrame(characteristic = c("Dataset Name", "Classifier Name"),
#'                        value = c("Asthma", "Different Means"))
#'              )
#'   #}
#' 
#' @export
setGeneric("runTests", function(measurements, ...)
           standardGeneric("runTests"))

setMethod("runTests", c("matrix"), # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  if(is.null(colnames(measurements)))
    stop("'measurements' matrix must have sample identifiers as its column names.")
  runTests(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("runTests", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, crossValParams = CrossValParams(), modellingParams = ModellingParams(),
                   characteristics = DataFrame(), verbose = 1)
{
  # Get out the classes if inside of data table.           
  if(is.null(rownames(measurements)))
    stop("'measurements' DataFrame must have sample identifiers as its row names.")
  
  if(any(is.na(measurements)))
    stop("Some data elements are missing and classifiers don't work with missing data. Consider imputation or filtering.")            
            
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  classes <- splitDataset[["classes"]]
  
  # Elements of list returned by runTest, in order.
  resultTypes <- c("ranked", "selected", "models", "testSet", "predictions", "tune")
  
  if("prevalidated" %in% names(modellingParams@trainParams))
  { # Check that all assays have a specification of how to classify them.
    prevalidate <- TRUE
    assaysWithoutParams <- setdiff(S4Vectors::mcols(measurements)[, "dataset"], c(names(modellingParams@trainParams), "clinical"))
    if(length(assaysWithoutParams) > 0)
      stop(paste("'trainParams' lacks a list for", paste(assaysWithoutParams, collapse = ", "), "assay."))
    return("To do.")
  } else {
    prevalidate <- FALSE
  }  # Re-implement pre-validation the right way later.
  
  featureInfo <- .summaryFeatures(measurements, prevalidate)
  allFeatures <- featureInfo[[1]]
  featureNames <- featureInfo[[2]]
  consideredFeatures <- featureInfo[[3]]
  # Create all partitions of training, testing and optionally validation sets.
  samplesSplits <- .samplesSplits(crossValParams, classes)
  splitsTestInfo <- .splitsTestInfo(crossValParams, samplesSplits)
  
  results <- bpmapply(function(trainingSamples, testSamples, setNumber)
  #results <- mapply(function(trainingSamples, testSamples, setNumber)
  {
    if(verbose >= 1 && setNumber %% 10 == 0)
      message("Processing sample set ", setNumber, '.')
    
    # crossValParams is needed at least for nested feature tuning.
    runTest(measurements, classes, trainingSamples, testSamples,
            crossValParams, modellingParams, characteristics, verbose,
            .iteration = setNumber)
  }, samplesSplits[["train"]], samplesSplits[["test"]], (1:length(samplesSplits[["train"]])),
  BPPARAM = crossValParams@parallelParams, SIMPLIFY = FALSE)
  #SIMPLIFY = FALSE)

  # Error checking and reporting.
  resultErrors <- sapply(results, function(result) is.character(result))
  if(sum(resultErrors) == length(results))
  {
      message("Error: All cross-validations had an error.")
      if(length(unique(unlist(results))) == 1)
        message("The common problem is: ", unlist(results)[[1]])
      return(results)
  } else if(sum(resultErrors) != 0) # Filter out cross-validations resulting in error.
  {
    warning(paste(sum(resultErrors),  "cross-validations, but not all, had an error and have been removed from the results."))
    results <- results[!resultErrors]
  }
  
  validationText <- .validationText(crossValParams)
  
  modParamsList <- list(modellingParams@transformParams, modellingParams@selectParams, modellingParams@trainParams, modellingParams@predictParams)
  autoCharacteristics <- lapply(modParamsList, function(stageParams) if(!is.null(stageParams)) stageParams@characteristics)
  autoCharacteristics <- do.call(rbind, autoCharacteristics)

  # Add extra settings which don't create varieties.
  # extras <- do.call(c, lapply(modParamsList, function(stageParams) if(!is.null(stageParams)) stageParams@otherParams))
  # if(length(extras) > 0)
  #   extras <- extras[sapply(extras, is.atomic)] # Store basic variables, not complex ones.
  # extrasDF <- DataFrame(characteristic = names(extras), value = unlist(extras))
  # characteristics <- rbind(characteristics, extrasDF)
  characteristics <- .filterCharacteristics(characteristics, autoCharacteristics)
  characteristics <- rbind(characteristics,
                             S4Vectors::DataFrame(characteristic = "Cross-validation", value = validationText))

  if(is.factor(results[[1]][["predictions"]]) | is.numeric(results[[1]][["predictions"]]))
    predictionsTable <- data.frame(sample = unlist(lapply(results, "[[", "testSet")), splitsTestInfo, class = unlist(lapply(results, "[[", "predictions")), check.names = FALSE)
  else # data frame
    predictionsTable <- data.frame(sample = unlist(lapply(results, "[[", "testSet")), splitsTestInfo, do.call(rbind, lapply(results, "[[", "predictions")), check.names = FALSE)
  rownames(predictionsTable) <- NULL
  tuneList <- lapply(results, "[[", "tune")
  if(length(unlist(tuneList)) == 0)
    tuneList <- NULL
  
  ClassifyResult(characteristics, rownames(measurements), allFeatures,
                 lapply(results, "[[", "ranked"), lapply(results, "[[", "selected"),
                 lapply(results, "[[", "models"), tuneList, predictionsTable, classes)
})

setMethod("runTests", c("MultiAssayExperiment"),
          function(measurements, targets = names(measurements), classes, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes, restrict = NULL)
  runTests(tablesAndClasses[["dataTable"]], tablesAndClasses[["classes"]], ...)            
})