setGeneric("easyHardClassifierTrain", function(measurements, ...)
standardGeneric("easyHardClassifierTrain"))

setMethod("easyHardClassifierTrain", "MultiAssayExperiment",
function(measurements, easyDatasetID = "clinical", hardDatasetID = names(measurements)[1],
         featureSets = NULL, metaFeatures = NULL, minimumOverlapPercent = 80,
         datasetName = NULL, classificationName = "Easy-Hard Classifier",
         easyClassifierParams = list(minCardinality = 10, minPurity = 0.9),
         hardClassifierParams = list(SelectParams(), TrainParams(), PredictParams()), 
         verbose = 3)
{
  if(easyDatasetID == "clinical")
  {
    easyDataset <- MultiAssayExperiment::colData(measurements) # Will be DataFrame
    easyDataset <- easyDataset[!is.na(easyDataset[, "class"]), ]
    easyDataset <- easyDataset[, -match("class", colnames(easyDataset))]
  } else if(easyDatasetID %in% names(measurements))
  {
    easyDataset <- measurements[, , easyDatasetID][[1]] # Get the underlying data container e.g. matrix.
    if(is.matrix(easyDataset))
      easyDataset <- t(easyDataset) # Make the variables be in columns.
  } else {
    stop("'easyDatasetID' is not \"clinical\" nor the name of any assay in 'measurements'.")
  }
  
  hardDataset <- measurements[, , hardDatasetID][[1]] # Get the underlying data container e.g. matrix.
  hardDataset <- S4Vectors::DataFrame(t(hardDataset), check.names = FALSE) # Variables as columns.
  commonSamples <- intersect(rownames(easyDataset), rownames(hardDataset))
  easyDataset <- easyDataset[commonSamples, ]
  hardDataset <- hardDataset[commonSamples, ]
  
  datasetIDs <- setNames(c(easyDatasetID, hardDatasetID), c("easy", "hard"))
  
  if(!requireNamespace("zoo", quietly = TRUE))
    stop("The package 'zoo' could not be found. Please install it.")
  
  if(verbose == 3)
    message("Fitting easy classifier to data.")
  
  classes <- MultiAssayExperiment::colData(measurements)[, "class"]
  names(classes) <- rownames(MultiAssayExperiment::colData(measurements))
  classes <- classes[match(rownames(easyDataset), names(classes))]
  easyDataset <- easyDataset[!is.na(classes), ]
  classes <- na.omit(classes)

  predictiveRules <- lapply(seq_along(easyDataset), function(featureIndex)
  { # Use a list rather than a table so that values can retain their type, such as numeric and categorical.
    measurement <- easyDataset[, featureIndex]
    if(is.numeric(measurement))
    {
      ordering <- order(measurement)
      midpoints <- zoo::rollmean(unique(na.omit(measurement[ordering])), 2) # Better performance if many repeated values.
      
      lowerRules <- lapply(midpoints, function(midpoint)
      {
        lowerClasses <- na.omit(classes[measurement < midpoint])
        lowerClassesCounts <- table(lowerClasses)
        lowerClassesProportions <- lowerClassesCounts / sum(lowerClassesCounts)
        isLarge <- length(lowerClasses) >= easyClassifierParams[["minCardinality"]]
        predictClass <- NULL
        if(any(lowerClassesProportions >= easyClassifierParams[["minPurity"]]))
          predictClass <- names(lowerClassesCounts)[which.max(lowerClassesCounts)]
        if(isLarge && !is.null(predictClass))
        {
          list(feature = colnames(easyDataset)[featureIndex],
               relation = '<', value = midpoint, predict = predictClass)
        }
      })
      whichRules <- which(sapply(lowerRules, length) > 0)
      if(any(whichRules)) # There is a worthwhile split.
        lowerRules <- lowerRules[whichRules[length(whichRules)]] # Get the rule with the biggest number of samples.
      else
        lowerRules <- NULL
      
      higherRules <- lapply(midpoints, function(midpoint)
      {
        higherClasses <- na.omit(classes[measurement > midpoint])
        higherClassesCounts <- table(higherClasses)
        higherClassesProportions <- higherClassesCounts / sum(higherClassesCounts)
        isLarge <- length(higherClasses) >= easyClassifierParams[["minCardinality"]]
        predictClass <- NULL
        if(any(higherClassesProportions >= easyClassifierParams[["minPurity"]]))
          predictClass <- names(higherClassesProportions)[which.max(higherClassesProportions)]
        if(isLarge && !is.null(predictClass))
        {
          list(feature = colnames(easyDataset)[featureIndex], relation = '>', value = midpoint, predict = predictClass)
        }
      })
      whichRules <- which(sapply(higherRules, length) > 0)
      if(any(whichRules)) # There is a worthwhile split.
        higherRules <- higherRules[whichRules[1]] # Get the rule with the biggest number of samples.
      else
        higherRules <- NULL

      c(lowerRules, higherRules)
    } else { # The variable is categorical.
      classesCounts <- table(measurement, classes)
      whichPureGroups <- which(classesCounts >= easyClassifierParams[["minCardinality"]] & t(apply(classesCounts, 1, function(measurementClasses) measurementClasses / sum(measurementClasses) > easyClassifierParams[["minPurity"]])), arr.ind = TRUE)
      if(nrow(whichPureGroups) > 0)
      {
        rownames(whichPureGroups) <- NULL
        apply(whichPureGroups, 1, function(valueClassPair)
        {
          list(feature = colnames(easyDataset)[featureIndex],
               relation = "==", value = rownames(classesCounts)[valueClassPair["measurement"]],
               predict = colnames(classesCounts)[valueClassPair["classes"]])          
        })
      }
    }
  })
  predictiveRules <- unlist(predictiveRules, recursive = FALSE)
  if(length(predictiveRules) > 0)
  {
    predictionsAndSamples <- .getEasyPredictions(easyDataset, predictiveRules)
    samplesHard <- predictionsAndSamples[[2]]
  } else { # No clinical variable is useful and all samples should be included for the hard classifier.
    samplesHard <- rownames(easyDataset)
  }
  
  if(!hardDatasetID %in% names(measurements))
  {
    stop("'hardDatasetID' is not the name of any assay in 'measurements'.")
  } else {
    if(length(samplesHard) == 0)
    {
      return(EasyHardClassifier(predictiveRules, NULL, datasetIDs))
    }

    hardDataset <- hardDataset[samplesHard, ]
    hardClasses <- classes[samplesHard]
    
    samplesHardClasses <- table(hardClasses)
    if(max(samplesHardClasses) >= sum(samplesHardClasses) - 1) # All samples or all samples except one belong to a particular class. Predict that class.
    {
      majorClass <- names(samplesHardClasses)[which.max(samplesHardClasses)]
      if(verbose == 3)
        message("Predicting ", length(samplesHard), " remaining samples as ", majorClass, '.')
      EasyHardClassifier(predictiveRules, list(selected = NULL, model = majorClass), datasetIDs)
    } else { # Train a classifier on the hard to classify by rules samples.
      if(verbose == 3)
        message("Fitting hard classifier to ", length(samplesHard), " samples not easily classified by easy classifier.")
      selectionAndModel <- runTest(hardDataset, hardClasses, "none", featureSets, metaFeatures, 80, datasetName, classificationName, names(hardClasses), names(hardClasses), hardClassifierParams, .iteration = "internal")
      if(is.character(selectionAndModel)) return(selectionAndModel) # An error occurred. One such case is when only one sample or no samples are left in a particular class.
      initialClass <- class(selectionAndModel[["models"]]) 
      if(!"list" %in% initialClass)
      {
        selectionAndModel[["selected"]] <- list(selectionAndModel[["selected"]])
        selectionAndModel[["models"]] <- list(selectionAndModel[["models"]])
      }
      trainedModels <- mapply(function(selected, model)
      {
        hardClassifier <- list(selected = selected, model = model)
        EasyHardClassifier(predictiveRules, hardClassifier, datasetIDs)
      }, selectionAndModel[["selected"]], selectionAndModel[["models"]], SIMPLIFY = FALSE)
      names(trainedModels) <- names(selectionAndModel[["models"]])
      
      if(initialClass != "list")
        trainedModels <- trainedModels[[1]]
      trainedModels
    }
  }
})

setGeneric("easyHardClassifierPredict", function(model, test, ...)
           standardGeneric("easyHardClassifierPredict"))

setMethod("easyHardClassifierPredict", c("EasyHardClassifier", "MultiAssayExperiment"), function(model, test, predictParams, verbose = 3)
{
  easyDatasetID <- model@datasetIDs["easy"]
  hardDatasetID <- model@datasetIDs["hard"]
  if(easyDatasetID == "clinical")
  {
    easyDataset <- MultiAssayExperiment::colData(test) # Will be DataFrame
  } else if(easyDataset %in% names(test))
  {
    easyDataset <- measurements[, , easyDataset][[1]] # Get the underlying data container e.g. matrix.
    if(is.matrix(easyDataset))
      easyDataset <- t(easyDataset) # Make the variables be in columns.
  }

  if(!is.null(model@easyClassifier))
  {
    predictionsAndSamples <- .getEasyPredictions(easyDataset, model@easyClassifier)
    samplesHard <- predictionsAndSamples[[2]]
  } else { # No clinical variable is useful and all samples should be included for the hard classifier.
    samplesHard <- rownames(easyDataset)
  }

  if(length(samplesHard) > 0)
  {
    if(!is.character(model@hardClassifier[["model"]]))
    {
      hardDataset <- test[, samplesHard, hardDatasetID][[1]] # Get the underlying data container e.g. matrix.
      hardDataset <- S4Vectors::DataFrame(t(hardDataset), check.names = FALSE) # Variables as columns.
      hardPredictions <- .doTest(model@hardClassifier[["model"]], hardDataset, samplesHard, predictParams, verbose)
    } else {
      hardPredictions <- rep(model@hardClassifier[["model"]], length(samplesHard)) # In training, samples belonging to one class left for hard classification.
    }
    
    predictionsClass <- class(hardPredictions)
    if(predictionsClass != "list")
      hardPredictions <- list(hardPredictions)
    
    allSamplesClasses <- lapply(hardPredictions, function(varietyPredictions)
    {
      allPredictions <- rep(NA, length(varietyPredictions)) # In case no rules are found for easy data set.
      if(exists("predictionsAndSamples")) # Rules exist for easy data set.
        allPredictions <- predictionsAndSamples[[1]]
      allPredictions[is.na(allPredictions)] <- as.character(varietyPredictions)
      allPredictions <- factor(allPredictions, levels = levels(MultiAssayExperiment::colData(test)[, "class"]))
      allPredictions    
    })

    if(predictionsClass != "list")
      allSamplesClasses <- unlist(allSamplesClasses, recursive = FALSE)
    allSamplesClasses    
  } else { # All samples predicted with easy classifier. Return those predictions.
    factor(predictionsAndSamples[[1]], levels(MultiAssayExperiment::colData(test)[, "class"]))
  }
})

.getEasyPredictions <- function(easyDataset, predictiveRules)
{
  # Now, determine the samples which are left to be classified by the hard classifier.
  allPredictions <- list()
  for(index in seq_along(predictiveRules))
  {
    value <- predictiveRules[[index]][["value"]]
    if(is.character(value))
      value <- paste('"', value, '"', sep = '')
    predictions <- rep(NA, nrow(easyDataset))
    isRuleTrue <- eval(parse(text = paste("easyDataset[, predictiveRules[[index]][[\"feature\"]]]", ' ', predictiveRules[[index]][["relation"]], ' ', value, sep = '')))
    predictions[isRuleTrue] <- predictiveRules[[index]][["predict"]]
    allPredictions[[index]] <- predictions
  }

  allPredictions <- do.call(cbind, allPredictions)
  samplesPredictionsCounts <- list()
  for(sampleIndex in 1:nrow(allPredictions)) # Always have some non-zero length for samplesPredictionsCounts.
  {
    samplesPredictionsCounts[[sampleIndex]] <- table(allPredictions[sampleIndex, ])
  }
  samplesPredictedClasses <- sapply(samplesPredictionsCounts, function(predictionsCounts)
  {
    if(length(predictionsCounts) == 0)
    {
      NULL
    } else if(length(predictionsCounts) == 1)
    {
      names(predictionsCounts)
    } else { # Two or more classes. Pick the most popular unless there is a tie.
      predictedClass <- names(predictionsCounts[which.max(predictionsCounts)])
      if(length(predictedClass) == 1)
        predictedClass
      else
        NULL
    }
  })
  samplesPredictedClasses[sapply(samplesPredictedClasses, is.null)] <- NA
  samplesHard <- rownames(easyDataset)[sapply(samplesPredictedClasses, is.na)]
  list(unlist(samplesPredictedClasses), samplesHard)
}