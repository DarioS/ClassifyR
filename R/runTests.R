# Reproducibly Run Various Kinds of Cross-Validation

setGeneric("runTests", function(measurements, ...) standardGeneric("runTests"))

setMethod("runTests", c("matrix"), function(measurements, outcome, ...) # Matrix of numeric measurements.
{
  if(is.null(rownames(measurements)))
    stop("'measurements' matrix must have sample identifiers as its row names.")
  runTests(S4Vectors::DataFrame(measurements, check.names = FALSE), outcome, ...)
})

setMethod("runTests", "DataFrame", function(measurements, outcome, crossValParams = CrossValParams(), modellingParams = ModellingParams(),
           characteristics = S4Vectors::DataFrame(), verbose = 1)
{
  # Get out the outcome if inside of data table.           
  if(is.null(rownames(measurements)))
    stop("'measurements' DataFrame must have sample identifiers as its row names.")
  
  if(any(is.na(measurements)))
    stop("Some data elements are missing and classifiers don't work with missing data. Consider imputation or filtering.")            

  originalFeatures <- colnames(measurements)
  if("feature" %in% colnames(mcols(measurements))) originalFeatures <- mcols(measurements)[, c("assay", "feature")]                 
  splitDataset <- prepareData(measurements, outcome)
  measurements <- splitDataset[["measurements"]]
  outcome <- splitDataset[["outcome"]]
  
  if(!is.null(modellingParams@selectParams) && max(modellingParams@selectParams@tuneParams[["nFeatures"]]) > ncol(measurements))
  {
      warning("Attempting to evaluate more features for feature selection than in
input data. Autmomatically reducing to smaller number.")
      modellingParams@selectParams@tuneParams[["nFeatures"]] <- 1:min(10, ncol(measurements))
  }
  
  # Element names of the list returned by runTest, in order.
  resultTypes <- c("ranked", "selected", "models", "testSet", "predictions", "tune", "importance")
  
  # Create all partitions of training and testing sets.
  samplesSplits <- .samplesSplits(crossValParams, outcome)
  splitsTestInfo <- .splitsTestInfo(crossValParams, samplesSplits)
  
  # Necessary hack for parallel processing on Windows.
  modellingParams <- modellingParams
  crossValParams <- crossValParams
  characteristics <- characteristics
  verbose <- verbose
  # Make them all local variables, so they are passed to workers.
  
  results <- bpmapply(function(trainingSamples, testSamples, setNumber)
  #results <- mapply(function(trainingSamples, testSamples, setNumber)
  {
    if(verbose >= 1 && setNumber %% 10 == 0)
      message("Processing sample set ", setNumber, '.')
    
    # crossValParams is needed at least for nested feature tuning.
    runTest(measurements[trainingSamples, , drop = FALSE], outcome[trainingSamples],
            measurements[testSamples, , drop = FALSE], outcome[testSamples],
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
  autoCharacteristics <- lapply(modParamsList, function(stageParams) if(!is.null(stageParams) && !is(stageParams, "PredictParams")) stageParams@characteristics)
  autoCharacteristics <- do.call(rbind, autoCharacteristics)

  # Add extra settings.
  extras <- unlist(lapply(modParamsList, function(stageParams) if(!is.null(stageParams)) stageParams@otherParams), recursive = FALSE)
  if(length(extras) > 0)
    extras <- extras[sapply(extras, is.atomic)] # Store basic variables, not complex ones.
  extrasDF <- S4Vectors::DataFrame(characteristic = names(extras), value = unname(unlist(extras)))
  characteristics <- rbind(characteristics, extrasDF)
  characteristics <- .filterCharacteristics(characteristics, autoCharacteristics)
  characteristics <- rbind(characteristics,
                             S4Vectors::DataFrame(characteristic = "Cross-validation", value = validationText))

  if(!is.data.frame(results[[1]][["predictions"]]))
  {
    if(is.numeric(results[[1]][["predictions"]])) # Survival task.
        predictsColumnName <- "risk"
    else # Classification task. A factor.
        predictsColumnName <- "class"
    predictionsTable <- S4Vectors::DataFrame(sample = unlist(lapply(results, "[[", "testSet")), splitsTestInfo, unlist(lapply(results, "[[", "predictions")), check.names = FALSE)
    colnames(predictionsTable)[ncol(predictionsTable)] <- predictsColumnName
  } else { # data frame
    predictionsTable <- S4Vectors::DataFrame(sample = unlist(lapply(results, "[[", "testSet")), splitsTestInfo, do.call(rbind, lapply(results, "[[", "predictions")), check.names = FALSE)
  }
  rownames(predictionsTable) <- NULL
  tuneList <- lapply(results, "[[", "tune")
  if(length(unlist(tuneList)) == 0)
    tuneList <- NULL
  importance <- NULL
  if(!is.null(results[[1]][["importance"]]))
    importance <- do.call(rbind, lapply(results, "[[", "importance"))
  
  ClassifyResult(characteristics, rownames(measurements), originalFeatures,
                 lapply(results, "[[", "ranked"), lapply(results, "[[", "selected"),
                 lapply(results, "[[", "models"), tuneList, predictionsTable, outcome, importance, modellingParams)
})

setMethod("runTests", c("MultiAssayExperiment"),
          function(measurements, targets = names(measurements), outcomeColumns, ...)
{
  omicsTargets <- setdiff(targets, "clinical")
  if(length(omicsTargets) > 0)
  {
    if(any(anyReplicated(measurements[, , omicsTargets])))
      stop("Data set contains replicates. Please provide remove or average replicate observations and try again.")
  }
  
  tablesAndOutcome <- .MAEtoWideTable(measurements, targets, outcomeColumns, restrict = NULL)
  runTests(tablesAndOutcome[["dataTable"]], tablesAndOutcome[["outcome"]], ...)            
})
