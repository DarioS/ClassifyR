setGeneric("runTest", function(measurements, ...)
           standardGeneric("runTest"))

setMethod("runTest", "matrix", # Matrix of numeric measurements.
  function(measurements, classes, ...)
{
  if(is.null(colnames(measurements)))
    stop("'measurements' matrix must have sample identifiers as its column names.")    
  runTest(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("runTest", "DataFrame", # Clinical data or one of the other inputs, transformed.
function(measurements, classes, training, testing, crossValParams = CrossValParams(), # crossValParams used for tuning optimisation.
         modellingParams = ModellingParams(), characteristics = DataFrame(), verbose = 1, .iteration = NULL)
{
  if(is.null(.iteration)) # Not being called by runTests but by user. So, check the user input.
  {
    if(is.null(rownames(measurements)))
      stop("'measurements' DataFrame must have sample identifiers as its row names.")
    if(any(is.na(measurements)))
      stop("Some data elements are missing and classifiers don't work with missing data. Consider imputation or filtering.")                
    
    splitDataset <- .splitDataAndClasses(measurements, classes)
    # Rebalance the class sizes of the training samples by either downsampling or upsampling
    # or leave untouched if balancing is none.
    rebalanced <- .rebalanceTrainingClasses(measurements, classes, training, testing, modellingParams@balancing)
    measurements <- rebalanced[["measurements"]]
    classes <- rebalanced[["classes"]]
    training <- rebalanced[["training"]]
    testing <- rebalanced[["testing"]]
    # Testing and set is not rebalanced, only the indices are shifted because the number
    # of rows of measurements changes, so the row numbers are updated, but the individual samples in testing
    # do not change.
  }
  
  testingSamplesIDs <- rownames(measurements)[testing]
  # All input features.
  if(!is.null(S4Vectors::mcols(measurements)))
    allFeatures <- S4Vectors::mcols(measurements)
  else
    allFeatures <- colnames(measurements)
  
  if(!is.null(modellingParams@transformParams))
  {
    if(length(modellingParams@transformParams@intermediate) != 0)
      modellingParams@transformParams <- .addIntermediates(modellingParams@transformParams)
    
    measurements <- tryCatch(.doTransform(measurements, training, modellingParams@transformParams, verbose), error = function(error) error[["message"]])
    if(is.character(measurements)) return(measurements) # An error occurred.
  }

  rankedFeatures <- NULL
  selectedFeatures <- NULL
  tuneDetailsSelect <- NULL
  if(!is.null(modellingParams@selectParams))
  {
    if(length(modellingParams@selectParams@intermediate) != 0)
      modellingParams@selectParams <- .addIntermediates(modellingParams@selectParams)
    
    topFeatures <- tryCatch(.doSelection(measurements, classes, training, modellingParams, verbose = verbose, crossValParams = crossValParams),
                            error = function(error) error[["message"]]) 

    if(is.character(topFeatures)) return(topFeatures) # An error occurred.
    rankedFeatures <- topFeatures[[1]] # Extract for result object.
    selectedFeatures <- topFeatures[[2]] # Extract for subsetting.
    tuneDetailsSelect <- topFeatures[[3]]
  
    if(modellingParams@selectParams@subsetToSelections == TRUE)
    { # Subset the the data table to only the selected features.
      if(is.null(S4Vectors::mcols(measurements)))
      { # Input was ordinary matrix or DataFrame.
          measurements <- measurements[, selectedFeatures, drop = FALSE]
      } else { # Input was MultiAssayExperiment. # Match the selected features to the data frame columns
          selectedIDs <-  do.call(paste, selectedFeatures)
          featuresIDs <- do.call(paste, S4Vectors::mcols(measurements)[, c("dataset", "feature")])
          selectedColumns <- match(selectedIDs, featuresIDs)
          measurements <- measurements[, selectedColumns, drop = FALSE]
      }
    }
  } 
  
  # Training stage.
  if(length(modellingParams@trainParams@intermediate) > 0)
    modellingParams@trainParams <- .addIntermediates(modellingParams@trainParams)
  
  trained <- tryCatch(.doTrain(measurements, classes, training, testing, modellingParams, verbose),
                      error = function(error) error[["message"]])

  if(is.character(trained)) return(trained) # An error occurred.
    
  tuneDetailsTrain <- trained[[2]] # Second element is tuning results.
  if(!is.null(modellingParams@trainParams@getFeatures)) # Features chosen inside classifier.
  {
    extrasList <- list()
    extras <- .methodFormals(modellingParams@trainParams@getFeatures)[-1]
    if(length(extras) > 0)
      extrasList <- mget(names(extras))
    
    featureInfo <- do.call(modellingParams@trainParams@getFeatures, c(trained[[1]], extrasList))
    rankedFeatures <- featureInfo[[1]]
    selectedFeatures <- featureInfo[[2]]
  }
  
  if(!is.null(modellingParams@predictParams))
  {
    if(length(modellingParams@predictParams@intermediate) != 0)
      modellingParams@predictParams <- .addIntermediates(modellingParams@predictParams)
                             
    predictedClasses <- tryCatch(.doTest(trained[["model"]], measurements, testing, modellingParams@predictParams, verbose),
                                error = function(error) error[["message"]]
                                )
    if(is.character(predictedClasses)) # An error occurred.
      return(predictedClasses) # Return early.
    
  } else {
    predictedClasses <- trained[[1]]
  }
  
  if(is.null(modellingParams@predictParams)) models <- NULL else models <- trained[[1]] # One function for training and testing. Typically, the models aren't returned to the user, such as Poisson LDA implemented by PoiClaClu.
  if(!is.null(tuneDetailsSelect)) tuneDetails <- tuneDetailsSelect else tuneDetails <- tuneDetailsTrain
  if(!is.null(.iteration)) # This function was called by runTests.
  {
    list(ranked = rankedFeatures, selected = selectedFeatures, models = models, testSet = testingSamplesIDs, predictions = predictedClasses, tune = tuneDetails)
  } else { # runTest is being used directly, rather than from runTests. Create a ClassifyResult object.
    # Only one training, so only one tuning choice, which can be summarised in characteristics.
    modParamsList <- list(modellingParams@transformParams, modellingParams@selectParams, modellingParams@trainParams, modellingParams@predictParams)
    if(!is.null(tuneDetails)) characteristics <- rbind(characteristics, data.frame(characteristic = colnames(tuneDetails),
                                                                                   value = unlist(tuneDetails)))
    autoCharacteristics <- do.call(rbind, lapply(modParamsList, function(stageParams) if(!is.null(stageParams)) stageParams@characteristics))
    characteristics <- .filterCharacteristics(characteristics, autoCharacteristics)
    characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Cross-validation", value = "Independent Set"))

    extras <- lapply(modParamsList, function(stageParams) if(!is.null(stageParams)) stageParams@otherParams)
    extrasDF <- DataFrame(characteristic = names(extras), value = unlist(extras))
    characteristics <- rbind(characteristics, extrasDF)
    
    ClassifyResult(characteristics, rownames(measurements), allFeatures, list(rankedFeatures), list(selectedFeatures),
                   list(models), tuneDetails, data.frame(sample = testing, class = predictedClasses), classes)
  }  
})

setMethod("runTest", c("MultiAssayExperiment"),
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, restrict = NULL)
  runTest(tablesAndClasses[["dataTable"]], tablesAndClasses[["classes"]], ...)            
})