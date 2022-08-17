# Perform a Single Classification

setGeneric("runTest", function(measurementsTrain, ...)
           standardGeneric("runTest"))

setMethod("runTest", "matrix", # Matrix of numeric measurements.
  function(measurementsTrain, outcomeTrain, measurementsTest, outcomeTest, ...)
{
  runTest(measurementsTrain = S4Vectors::DataFrame(measurementsTrain, check.names = FALSE),
          outcomeTrain = outcomeTrain,
          measurementsTest = S4Vectors::DataFrame(measurementsTest, check.names = FALSE),
          outcomeTest = outcomeTest,
          ...)
})

setMethod("runTest", "DataFrame", # Sample information data or one of the other inputs, transformed.
function(measurementsTrain, outcomeTrain, measurementsTest, outcomeTest,
         crossValParams = CrossValParams(), # crossValParams might be used for tuning optimisation.
         modellingParams = ModellingParams(), characteristics = S4Vectors::DataFrame(), verbose = 1, .iteration = NULL)
{
  if(is.null(.iteration)) # Not being called by runTests but by user. So, check the user input.
  {
    if(is.null(rownames(measurementsTrain)))
      stop("'measurementsTrain' DataFrame must have sample identifiers as its row names.")
    if(any(is.na(measurementsTrain)))
      stop("Some data elements are missing and classifiers don't work with missing data. Consider imputation or filtering.")                
    
    splitDatasetTrain <- prepareData(measurementsTrain, outcomeTrain)
      
    # Rebalance the class sizes of the training samples by either downsampling or upsampling
    # or leave untouched if balancing is none.
    if(!is(outcomeTrain, "Surv"))
    {
      rebalancedTrain <- .rebalanceTrainingClasses(splitDatasetTrain[["measurements"]], splitDatasetTrain[["outcome"]], modellingParams@balancing)
      measurementsTrain <- rebalancedTrain[["measurementsTrain"]]
      outcomeTrain <- rebalancedTrain[["classesTrain"]]
    }
  }
  
  if("feature" %in% colnames(mcols(measurementsTrain))) originalFeatures <- mcols(measurementsTrain)[, na.omit(match(c("assay", "feature"), colnames(mcols(measurementsTrain))))]    
  else originalFeatures <- colnames(measurementsTrain)
    
  if(!is.null(modellingParams@selectParams) && max(modellingParams@selectParams@tuneParams[["nFeatures"]]) > ncol(measurementsTrain))
  {
    warning("Attempting to evaluate more features for feature selection than in
input data. Autmomatically reducing to smaller number.")
    modellingParams@selectParams@tuneParams[["nFeatures"]] <- 1:min(10, ncol(measurementsTrain))
  }
  
  if(!is.null(crossValParams) && !is.null(crossValParams@adaptiveResamplingDelta))
  { # Iteratively resample training samples until their class probability or risk changes, on average, less than delta.
    delta <- crossValParams@adaptiveResamplingDelta
    crossValParams@adaptiveResamplingDelta <- NULL
    scoresPrevious <- rep(1, nrow(measurementsTrain))
    repeat{
      newSamples <- sample(nrow(measurementsTrain), replace = TRUE, prob = scoresPrevious)
      measurementsTrainResampled <- measurementsTrain[newSamples, ]
      outcomeResampled <- outcomeTrain[newSamples]
      ASpredictions <- runTest(measurementsTrainResampled, outcomeResampled,
              measurementsTrain, outcomeTrain, crossValParams, modellingParams,
              .iteration = "internal")[["predictions"]]
      if(is.factor(outcomeResampled))
          scoresNew <- mapply(function(rowIndex, class) ASpredictions[rowIndex, class], 1:nrow(ASpredictions), as.character(outcomeTrain))
      else
          scoresNew <- ASpredictions[, "risk"]

      if(mean(abs(scoresNew - scoresPrevious)) < delta)
      {
        measurementsTrain <- measurementsTrainResampled
        break;
      }
      scoresPrevious <- scoresNew
    }
  }
  
  if(!is.null(modellingParams@transformParams))
  {
    if(length(modellingParams@transformParams@intermediate) != 0)
      modellingParams@transformParams <- .addIntermediates(modellingParams@transformParams)
    
    measurementsTransformedList <- tryCatch(.doTransform(measurementsTrain, measurementsTest, modellingParams@transformParams, verbose), error = function(error) error[["message"]])
    if(is.character(measurementsTransformedList[[1]])) return(measurementsTransformedList[[1]]) # An error occurred.
  }

  rankedFeaturesIndices <- NULL
  selectedFeaturesIndices <- NULL
  tuneDetailsSelect <- NULL

  if(!is.null(modellingParams@selectParams))
  {
    if(length(modellingParams@selectParams@intermediate) != 0)
      modellingParams@selectParams <- .addIntermediates(modellingParams@selectParams)
 
    topFeatures <- tryCatch(.doSelection(measurementsTrain, outcomeTrain, crossValParams, modellingParams, verbose),
                            error = function(error) error[["message"]]) 
    if(is.character(topFeatures)) return(topFeatures) # An error occurred.
    
    rankedFeaturesIndices <- topFeatures[[1]] # Extract for result object.
    selectedFeaturesIndices <- topFeatures[[2]] # Extract for subsetting.
    tuneDetailsSelect <- topFeatures[[3]]

    if(modellingParams@selectParams@subsetToSelections == TRUE)
    {
      measurementsTrain <- measurementsTrain[, selectedFeaturesIndices, drop = FALSE]
      measurementsTest <- measurementsTest[, selectedFeaturesIndices, drop = FALSE]
    }
  }

  # Training stage.
  if(length(modellingParams@trainParams@intermediate) > 0)
    modellingParams@trainParams <- .addIntermediates(modellingParams@trainParams)
  if(!is.null(tuneDetailsSelect))
  {
    avoidTune <- match(colnames(tuneDetailsSelect), names(modellingParams@trainParams@tuneParams))
    if(any(!is.na(avoidTune)))
    {
      modellingParams@trainParams@otherParams <- c(modellingParams@trainParams@otherParams, tuneDetailsSelect[!is.na(avoidTune)])
      modellingParams@trainParams@tuneParams <- modellingParams@trainParams@tuneParams[-na.omit(avoidTune)]
      if(length(modellingParams@trainParams@tuneParams) == 0) modellingParams@trainParams@tuneParams <- NULL
    }
  }
  
  # Some classifiers have one function for training and testing, so that's why test data is also passed in.
  trained <- tryCatch(.doTrain(measurementsTrain, outcomeTrain, measurementsTest, outcomeTest, modellingParams, verbose),
                      error = function(error) error[["message"]])
  if(is.character(trained)) return(trained) # An error occurred.
  
  tuneDetailsTrain <- trained[[2]] # Second element is tuning results.
  
  if(!is.null(modellingParams@trainParams@getFeatures)) # Features chosen inside classifier.
  {
    extrasList <- list()
    extras <- .methodFormals(modellingParams@trainParams@getFeatures, class(trained[[1]]))[-1]
    if(length(extras) > 0)
      extrasList <- mget(setdiff(names(extras), "..."))

    rankedChosenList <- do.call(modellingParams@trainParams@getFeatures, c(trained[1], extrasList))
    rankedFeaturesIndices <- rankedChosenList[[1]]
    selectedFeaturesIndices <- rankedChosenList[[2]]
  }
  
  if(!is.null(modellingParams@predictParams))
  {
    if(length(modellingParams@predictParams@intermediate) != 0)
      modellingParams@predictParams <- .addIntermediates(modellingParams@predictParams)
    
    predictedOutcome <- tryCatch(.doTest(trained[["model"]], measurementsTest, modellingParams@predictParams, verbose),
                                error = function(error) error[["message"]]
                                )

    if(is.character(predictedOutcome)) # An error occurred.
      return(predictedOutcome) # Return early.
    
  } else { # One function that does training and testing, so predictions were made earlier
           # by .doTrain, rather than this .doTest stage.
    predictedOutcome <- trained[[1]]
  }
  
  # Exclude one feature at a time, build model, predict test samples.
  importanceTable <- NULL
  if(is.numeric(.iteration) && modellingParams@doImportance == TRUE)
  {
    performanceMP <- modellingParams@selectParams@tuneParams[["performanceType"]]
    performanceType <- ifelse(!is.null(performanceMP), performanceMP, "Balanced Error")
    performancesWithoutEach <- sapply(selectedFeaturesIndices, function(selectedIndex)
    {
      measurementsTrainLess1 <- measurementsTrain[, -selectedIndex, drop = FALSE]
      measurementsTestLess1 <- measurementsTest[, -selectedIndex, drop = FALSE]
      modelWithoutOne <- tryCatch(.doTrain(measurementsTrainLess1, outcomeTrain, measurementsTestLess1, outcomeTest, modellingParams, verbose),
                                  error = function(error) error[["message"]])
      if(!is.null(modellingParams@predictParams))
      predictedOutcomeWithoutOne <- tryCatch(.doTest(modelWithoutOne[["model"]], measurementsTestLess1, modellingParams@predictParams, verbose),
                                              error = function(error) error[["message"]])
      else predictedOutcomeWithoutOne <- modelWithoutOne[["model"]]

      if(!is.null(ncol(predictedOutcomeWithoutOne)))
        predictedOutcomeWithoutOne <- predictedOutcomeWithoutOne[, na.omit(match(c("class", "risk"), colnames(predictedOutcomeWithoutOne)))]
      calcExternalPerformance(outcomeTest, predictedOutcomeWithoutOne, performanceType)
    })
    
    if(!is.null(ncol(predictedOutcome)))
        predictedOutcome <- predictedOutcome[, na.omit(match(c("class", "risk"), colnames(predictedOutcome)))]
    performanceChanges <- round(performancesWithoutEach - calcExternalPerformance(outcomeTest, predictedOutcome, performanceType), 2)
     
    if(is.null(S4Vectors::mcols(measurementsTrain)))
    {
      selectedFeatures <- colnames(measurementsTrain)[selectedFeaturesIndices]
    } else {
      featureColumns <- na.omit(match(c("assay", "feature"), colnames(S4Vectors::mcols(measurementsTrain))))
      selectedFeatures <- S4Vectors::mcols(measurementsTrain)[selectedFeaturesIndices, featureColumns]
    }
    importanceTable <- S4Vectors::DataFrame(selectedFeatures, performanceChanges)
    if(ncol(importanceTable) == 2) colnames(importanceTable)[1] <- "feature"
    colnames(importanceTable)[ncol(importanceTable)] <- paste("Change in", performanceType)
  }
  
  if(is.null(modellingParams@predictParams)) models <- NULL else models <- trained[[1]] # One function for training and testing. Typically, the models aren't returned to the user, such as Poisson LDA implemented by PoiClaClu.
  if(!is.null(tuneDetailsSelect)) tuneDetails <- tuneDetailsSelect else tuneDetails <- tuneDetailsTrain

  # Convert back into original, potentially unsafe feature identifiers unless it is a nested cross-validation.
  if(is.null(.iteration) || .iteration != "internal")
  {
    if(!is.null(rankedFeaturesIndices))
    {
      if(is.null(S4Vectors::mcols(measurementsTrain)))
      {
        rankedFeatures <- originalFeatures[rankedFeaturesIndices]
      } else {
        featureColumns <- na.omit(match(c("assay", "feature"), colnames(S4Vectors::mcols(measurementsTrain))))          
        rankedFeatures <- originalFeatures[rankedFeaturesIndices, featureColumns]
      }
    } else { rankedFeatures <- NULL}
    if(!is.null(selectedFeaturesIndices))
    {
      if(is.null(S4Vectors::mcols(measurementsTrain))){
        selectedFeatures <- originalFeatures[selectedFeaturesIndices]
      } else {
        featureColumns <- na.omit(match(c("assay", "feature"), colnames(S4Vectors::mcols(measurementsTrain))))  
        selectedFeatures <- originalFeatures[selectedFeaturesIndices, ]
      }
    } else { selectedFeatures <- NULL}
  } else { # Nested use in feature selection. No feature selection in inner execution, so ignore features. 
      rankedFeatures <- selectedFeatures <- NULL
  }
  
  if(!is.null(.iteration)) # This function was not called by the end user.
  {
    list(ranked = rankedFeatures, selected = selectedFeatures, models = models, testSet = rownames(measurementsTest), predictions = predictedOutcome, tune = tuneDetails, importance = importanceTable)
  } else { # runTest executed by the end user. Create a ClassifyResult object.
    # Only one training, so only one tuning choice, which can be summarised in characteristics.
    modParamsList <- list(modellingParams@transformParams, modellingParams@selectParams, modellingParams@trainParams, modellingParams@predictParams)
    if(!is.null(tuneDetails)) characteristics <- rbind(characteristics, data.frame(characteristic = colnames(tuneDetails),
                                                                                   value = unlist(tuneDetails)))
    autoCharacteristics <- do.call(rbind, lapply(modParamsList, function(stageParams) if(!is.null(stageParams) && !is(stageParams, "PredictParams")) stageParams@characteristics))
    characteristics <- .filterCharacteristics(characteristics, autoCharacteristics)
    characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Cross-validation", value = "Independent Set"))
    
    extras <- unlist(lapply(modParamsList, function(stageParams) if(!is.null(stageParams)) stageParams@otherParams), recursive = FALSE)
    if(length(extras) > 0)
        extras <- extras[sapply(extras, is.atomic)] # Store basic variables, not complex ones.
    extrasDF <- S4Vectors::DataFrame(characteristic = names(extras), value = unname(unlist(extras)))
    characteristics <- rbind(characteristics, extrasDF)
    
    allSamples <- c(rownames(measurementsTrain), rownames(measurementsTest))
    if(!is.null(ncol(outcomeTrain)))
    {
      allOutcome <- rbind(outcomeTrain, outcomeTest)
      rownames(allOutcome) <- allSamples
    } else { 
      allOutcome <- c(outcomeTrain, outcomeTest)
      names(allOutcome) <- allSamples
    }

    ClassifyResult(characteristics, allSamples, originalFeatures, list(rankedFeatures), list(selectedFeatures),
                   list(models), tuneDetails, S4Vectors::DataFrame(sample = rownames(measurementsTest), predictedOutcome, check.names = FALSE), allOutcome, importanceTable)
  }  
})

setMethod("runTest", c("MultiAssayExperiment"),
          function(measurementsTrain, measurementsTest, targets = names(measurements), outcomeColumns, ...)
{
  omicsTargets <- setdiff(targets, "clinical")
  if(length(omicsTargets) > 0)
  {
    if(any(anyReplicated(measurements[, , omicsTargets])))
      stop("Data set contains replicates. Please provide remove or average replicate observations and try again.")
  }
  
  tablesAndClassesTrain <- .MAEtoWideTable(measurementsTrain, targets, outcomeColumns, restrict = NULL)
  tablesAndClassesTest <- .MAEtoWideTable(measurementsTest, targets, outcomeColumns, restrict = NULL)
  runTest(tablesAndClassesTrain[["dataTable"]], tablesAndClassesTrain[["outcome"]],
          tablesAndClassesTest[["dataTable"]], tablesAndClassesTest[["outcome"]], ...)            
})
