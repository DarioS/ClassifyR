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
  extras <- do.call(c, lapply(modParamsList, function(stageParams) if(!is.null(stageParams)) stageParams@otherParams))
  if(length(extras) > 0)
    extras <- extras[sapply(extras, is.atomic)] # Store basic variables, not complex ones.
  extrasDF <- DataFrame(characteristic = names(extras), value = unlist(extras))
  characteristics <- rbind(characteristics, extrasDF)
  characteristics <- .filterCharacteristics(characteristics, autoCharacteristics)
  characteristics <- rbind(characteristics,
                             S4Vectors::DataFrame(characteristic = "Cross-validation", value = validationText))

  if(is.factor(results[[1]][["predictions"]]))
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