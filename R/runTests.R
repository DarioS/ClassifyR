setGeneric("runTests", function(measurements, ...)
           {standardGeneric("runTests")})

setMethod("runTests", c("matrix"), # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  if(is.null(colnames(measurements)))
    stop("'measurements' matrix must have sample identifiers as its column names.")
  runTests(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("runTests", c("DataFrame"), # Clinical data or one of the other inputs, transformed.
          function(measurements, classes,
                      featureSets = NULL, metaFeatures = NULL, minimumOverlapPercent = 80,
                      datasetName, classificationName,
                      validation = c("permute", "leaveOut", "fold"),
                      permutePartition = c("fold", "split"),
                      permutations = 100, percent = 25, folds = 5, leave = 2,
                      seed, parallelParams = bpparam(),
                      params = list(SelectParams(), TrainParams(), PredictParams()),
                      verbose = 1)
{
  if(is.null(rownames(measurements)))
    stop("'measurements' DataFrame must have sample identifiers as its row names.")
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  classes <- splitDataset[["classes"]]

  validation <- match.arg(validation)
  permutePartition <- match.arg(permutePartition)
  if(!missing(seed)) set.seed(seed)
  resultTypes <- c("ranked", "selected", "models", "testSet", "predictions", "tune")
  # Elements of list returned by runTest.

  stagesParamClasses <- sapply(params, class)
  if(match("TrainParams", stagesParamClasses) > match("PredictParams", stagesParamClasses))
    stop("\"testing\" must not be before \"training\" in 'params'.")
  
  selectParams <- params[[match("SelectParams", stagesParamClasses)]]
  predictParams <- params[[match("PredictParams", stagesParamClasses)]]
  
  if(!is.null(S4Vectors::mcols(measurements)))
  {
    allFeatures <- S4Vectors::mcols(measurements)
    featureNames <- S4Vectors::mcols(measurements)[, "feature"]
  } else {
    allFeatures <- colnames(measurements)
    featureNames <- colnames(measurements)
  }
  
  if(!is.null(featureSets)) # Feature sets provided and automatically filtered.
  {
    # Filter out the edges or the features from the sets which are not in measurements.
    featureSetsList <- featureSets@sets
    if(class(featureSetsList[[1]]) == "matrix")
    {
      setsSizes <- sapply(featureSetsList, nrow)
      edgesAll <- do.call(rbind, featureSetsList)
      networkNames <- rep(names(featureSetsList), setsSizes)
      edgesKeep <- edgesAll[, 1] %in% featureNames & edgesAll[, 2] %in% featureNames
      edgesFiltered <- edgesAll[edgesKeep, ]
      networkNamesFiltered <- networkNames[edgesKeep]
      setsRows <- split(1:nrow(edgesFiltered), factor(networkNamesFiltered, levels = names(featureSetsList)))
      featureSetsListFiltered <- lapply(setsRows, function(setRows) edgesFiltered[setRows, , drop = FALSE])
      setsSizesFiltered <- sapply(featureSetsListFiltered, nrow)
    } else { # A set of features without edges, such as a gene set.
      setsSizes <- sapply(setsNodes, length)
      nodesVector <- unlist(featureSetsList)
      setsVector <- rep(names(featureSetsList), setsSizes)
      keepNodes <- !is.na(match(nodesVector, featureNames))
      nodesVector <- nodesVector[keepNodes]
      setsVector <- setsVector[keepNodes]
      featureSetsListFiltered <- split(nodesVector, factor(setsVector, levels = names(featureSetsList)))
      setsSizesFiltered <- sapply(featureSetsListFiltered, length)
    }
      
    keepSets <- setsSizesFiltered / setsSizes * 100 >= minimumOverlapPercent
    featureSetsListFiltered <- featureSetsListFiltered[keepSets]
    featureSets <- FeatureSetCollection(featureSetsListFiltered)
    measurements <- measurements[, featureNames %in% unlist(featureSetsListFiltered)]
      
    if(verbose >= 1)
      message("After filtering features, ", length(featureSetsListFiltered), " out of ", length(featureSetsList), " sets remain.")
  }

  # Could refer to features or feature sets, depending on if a selection method utilising feature sets is used.
  if(!is.null(selectParams) && !is.null(featureSets))
    consideredFeatures <- length(featureSets@sets)
  else
    consideredFeatures <- ncol(measurements)
  
  if(validation == "permute")
  {
    if(permutePartition == "fold")
    {
      samplesFolds <- lapply(1:permutations, function(permutation)
                      {
                        classesFolds <- lapply(levels(classes), function(className)
                        {
                          whichSamples <- which(classes == className)
                          split(sample(whichSamples), rep(1:folds, length.out = length(whichSamples)))
                        })
                        allFolds <- lapply(1:folds, function(fold)
                        {
                          unlist(lapply(classesFolds, "[[", fold))
                        })
                      })
    } else { # Is split.
      samplesTrain <- (100 - percent) / 100 * table(classes)
      samplesFolds <- lapply(1:permutations, function(permutation)
                      {
                        trainSet <- unlist(mapply(function(className, number)
                        {
                          sample(which(classes == className), number)
                        }, levels(classes), samplesTrain))
                        testSet <- setdiff(1:length(classes), trainSet)
                        list(trainSet, testSet)
                      })
    }

    results <- bpmapply(function(sampleFolds, sampleNumber)
    {
      if(verbose >= 1 && sampleNumber %% 10 == 0)
        message("Processing sample set ", sampleNumber, '.')
      if(permutePartition == "fold")
      {
        lapply(1:folds, function(foldIndex)
        {
          runTest(measurements, classes, featureSets, metaFeatures, minimumOverlapPercent, training = unlist(sampleFolds[-foldIndex]),
                  testing = sampleFolds[[foldIndex]], params = params, verbose = verbose,
                  .iteration = c(sampleNumber, foldIndex))
        })
      } else { # Split mode.
        runTest(measurements, classes, featureSets, metaFeatures, minimumOverlapPercent, training = sampleFolds[[1]],
                testing = sampleFolds[[2]], params = params, verbose = verbose, .iteration = sampleNumber)
      }
    }, samplesFolds, as.list(1:permutations), BPPARAM = parallelParams, SIMPLIFY = FALSE)
  } else if(validation == "leaveOut") # leave k out.
  {
    testSamples <- as.data.frame(utils::combn(nrow(measurements), leave))
    trainingSamples <- lapply(testSamples, function(sample) setdiff(1:nrow(measurements), sample))
    results <- bpmapply(function(trainingSample, testSample, sampleNumber)
    {
      if(verbose >= 1 && sampleNumber %% 10 == 0)
        message("Processing sample set ", sampleNumber, '.')
      runTest(measurements, classes, featureSets, metaFeatures, minimumOverlapPercent, training = trainingSample, testing = testSample,
              params = params, verbose = verbose, .iteration = sampleNumber)
    }, trainingSamples, testSamples, (1:length(trainingSamples)),
    BPPARAM = parallelParams, SIMPLIFY = FALSE)
  } else { # Unresampled, ordinary k-fold cross-validation.
      classesFolds <- lapply(levels(classes), function(className)
      {
        whichSamples <- which(classes == className)
        split(sample(whichSamples), rep(1:folds, length.out = length(whichSamples)))
      })
      samplesFolds <- lapply(1:folds, function(fold)
      {
        unlist(lapply(classesFolds, "[[", fold))
      })
    
    if(verbose >= 1)
      message("Processing ", folds, "-fold cross-validation.")
    results <- bplapply(1:folds, function(foldIndex)
    {
      runTest(measurements, classes, featureSets, metaFeatures, minimumOverlapPercent, training = unlist(samplesFolds[-foldIndex]),
              testing = samplesFolds[[foldIndex]], params = params, verbose = verbose,
              .iteration = foldIndex)
    }, BPPARAM = parallelParams)
  }

  if(validation == "permute" && permutePartition == "split" || validation %in% c("leaveOut", "fold"))
  {
    resultErrors <- sapply(results, function(result) is.character(result))
    if(sum(resultErrors) == length(results))
    {
      message("Error: All cross-validations had an error.")
      if(length(unique(unlist(results))) == 1)
        message("The common problem is: ", unlist(results)[[1]])
      return(results)
    } else if(sum(resultErrors) != 0) # Filter out cross-validations resulting in error.
    {
      warning(paste(sum(resultErrors),  "cross-validations had an error and have been removed from the results."))
      results <- results[!resultErrors]
    }
  } else { # Result has nested lists, because of permutation reordering with folding.
    resultErrors <- lapply(results, function(resample) lapply(resample, is.character))
    if(sum(unlist(resultErrors)) == permutations * folds)
    {
      message("Error: All cross-validations had an error.")
      if(length(unique(unlist(results))) == 1)
        message("The common problem is: ", unlist(results)[[1]])
      return(results)
    } else if(sum(unlist(resultErrors)) != 0) # Filter out error cross-validations.
    {
      warning(paste(sum(unlist(resultErrors)),  "cross-validations had an error and have been removed from the results."))	    
      results <- results[sapply(results, function(resample) !all(sapply(resample, is.character)))]
      results <- lapply(results, function(resample) resample[!sapply(resample, is.character)])
    }
  }
  
  if(validation == "permute")
  {
    if(permutePartition == "fold")
    {
      if(class(results[[1]][[1]][["predictions"]]) == "list") 
      {
        multipleVarieties <- TRUE
        varietyNames <- names(results[[1]][[1]][["predictions"]])
      } else {
        multipleVarieties <- FALSE
        varietyNames <- "None"
      }
    } else {
      if(class(results[[1]][["predictions"]]) == "list")
      {
        multipleVarieties <- TRUE
        varietyNames <- names(results[[1]][["predictions"]])
      } else {
        multipleVarieties <- FALSE
        varietyNames <- "None"
      }
    }

    resultsByVariety <- lapply(varietyNames, function(varietyName)
    {
      resultsLists <- lapply(resultTypes, function(resultType)
      {
        resultList <- lapply(results, function(sample)
        {
          if(resultType == "testSet" || multipleVarieties == FALSE)
          {
            if(resultType %in% c("testSet", "predictions"))
            {
              if(permutePartition == "fold")
              {
                if(resultType == "testSet" || class(sample[[1]][["predictions"]]) != "data.frame")
                {
                  # First fold predictions are a not data.frame.
                  reshaped <- unlist(lapply(sample, function(fold) fold[[resultType]]))
                  if(resultType == "predictions")
                    reshaped <- factor(reshaped, levels = levels(classes))
                  reshaped
                } else { # Predictions is a data.frame with classes and scores.
                  resultTable <- do.call(rbind, lapply(sample, function(fold) fold[[resultType]]))
                  resultTable[, "class"] <- factor(resultTable[, "class"], levels = levels(classes))
                  resultTable
                }
              } else { # permutePartition is "split".
                information <- sample[[resultType]]
                if(resultType == "predictions")
                  information <- factor(information, levels = levels(classes))
                information
              }
            } else {
              if(permutePartition == "fold")
                lapply(sample, function(fold) fold[[resultType]])
              else
                sample[[resultType]]
            }
          } else { # There are multiple varieties.
            if(resultType == "predictions") # Predictions with varieties.
            {
              if(permutePartition == "fold")
              {
                prediction <- lapply(sample, function(fold) fold[[resultType]][[varietyName]])
                if(class(prediction[[1]]) == "data.frame")
                {
                  reshaped <- do.call(rbind, prediction)
                  reshaped[, "class"] <- factor(reshaped[, "class"], levels = levels(classes))
                  reshaped
                } else
                  factor(unlist(prediction), levels = levels(classes))
              } else { # permutePartition is "split".
                information <- sample[[resultType]][[varietyName]]
                if(resultType == "predictions")
                  information <- factor(information, levels = levels(classes))
                information                
              }
            } else { # Ranked, selected features, models or test samples, with varieties.
              if(permutePartition == "fold")
              {
                if(is.null(sample[[1]][[resultType]]))
                  NULL
                else
                  lapply(sample, function(fold) fold[[resultType]][[varietyName]])
              } else {
                if(is.null(sample[[resultType]]))
                  NULL
                else
                  sample[[resultType]][[varietyName]]
              }
            }
          }
        })
        if(is.null(unlist(resultList, use.names = FALSE))) list(NULL) else resultList # Collapse multiple NULLs into one NULL.
      })
      names(resultsLists) <- resultTypes
      resultsLists
    })
    names(resultsByVariety) <- varietyNames
  } else {# leave k out or ordinary, unresampled k-fold cross-validation.
    if(class(results[[1]][["predictions"]]) == "list")
    {
      multipleVarieties <- TRUE
      varietyNames <- names(results[[1]][["predictions"]])
    } else {
      multipleVarieties <- FALSE
      varietyNames <- "None"
    }
    
    resultsByVariety <- lapply(varietyNames, function(varietyName)
    {
      resultsLists <- lapply(resultTypes, function(resultType)
      {
        resultList <- lapply(results, function(sample) sample[[resultType]])
        if(resultType %in% c("testSet", "predictions") && multipleVarieties == FALSE && class(resultList[[1]]) != "data.frame")
        {
          reshaped <- unlist(resultList)
          if(resultType == "predictions")
            reshaped <- factor(reshaped, levels = levels(classes)) # Enforce that levels are in the same order as input by the user.
          reshaped
        } else if (resultType == "predictions" && multipleVarieties == FALSE) { # Predictions are in a data.frame.
          resultTable <- do.call(rbind, resultList)
          resultTable[, "class"] <- factor(resultTable[, "class"], levels = levels(classes))
          resultTable
        } else if(multipleVarieties == TRUE && resultType == "predictions" && class(resultList[[1]]) != "data.frame") {
          factor(unlist(lapply(resultList, "[[", varietyName)), levels = levels(classes))
        } else if(multipleVarieties == TRUE && resultType == "predictions") { # Predictions are factor or numeric.
          resultTable <- do.call(rbind, lapply(resultList, "[[", varietyName))
          resultTable[, "class"] <- factor(resultTable[, "class"], levels = levels(classes))
          resultTable
        } else if(multipleVarieties == TRUE && resultType %in% c("ranked", "selected", "tune", "models")) {
          returnValue <- lapply(resultList, "[[", varietyName)
          if(is.null(unlist(returnValue, use.names = FALSE))) NULL else returnValue
        } else { # For multipleVarieties being FALSE, and not the predicted classes.
          resultList
        }
      })
      names(resultsLists) <- resultTypes
      resultsLists
    })
    names(resultsByVariety) <- varietyNames
  }
  
  predictionTablesByVariety <- lapply(resultsByVariety, function(resultVariety)
  {
    if(validation == "permute")
    {
      lapply(1:length(results), function(resample)
      {
        sampleNames <- resultVariety[["testSet"]][[resample]]
        switch(class(resultVariety[["predictions"]][[resample]]),
               factor = data.frame(sample = sampleNames, class = factor(resultVariety[["predictions"]][[resample]], levels = levels(classes)), stringsAsFactors = FALSE, row.names = NULL),
               numeric = data.frame(sample = sampleNames, resultVariety[["predictions"]][[resample]], stringsAsFactors = FALSE, row.names = NULL),
               data.frame = data.frame(sample = sampleNames,
                                       class = factor(resultVariety[["predictions"]][[resample]][, sapply(resultVariety[["predictions"]][[resample]], class) == "factor"], levels = levels(classes)),
                                       resultVariety[["predictions"]][[resample]][, sapply(resultVariety[["predictions"]][[resample]], class) == "numeric"], stringsAsFactors = FALSE, row.names = NULL))
                
      })
    } else { # leave k out or ordinary, unresampled k-fold cross-validation.
      sampleNames <- resultVariety[["testSet"]]
      list(switch(class(resultVariety[["predictions"]]),
             factor = data.frame(sample = sampleNames, class = factor(resultVariety[["predictions"]]), stringsAsFactors = FALSE, row.names = NULL),
             numeric = data.frame(sample = sampleNames, resultVariety[["predictions"]], stringsAsFactors = FALSE, row.names = NULL),
             data.frame = data.frame(sample = sampleNames,
                                     class = factor(resultVariety[["predictions"]][, sapply(resultVariety[["predictions"]], class) == "factor"], levels = levels(classes)),
                                     resultVariety[["predictions"]][, sapply(resultVariety[["predictions"]], class) == "numeric"], stringsAsFactors = FALSE, row.names = NULL)))
    }
  })

  resultsByVariety <- mapply(function(resultVariety, predictionVariety) # Put formatted prediction tables into results list.
  {
    resultVariety[["predictions"]] <- predictionVariety
    resultVariety
  }, resultsByVariety, predictionTablesByVariety, SIMPLIFY = FALSE)
  
  if(validation == "permute")
  {
    if(permutePartition == "fold")
      validationInfo <- list("permuteFold", permutations, folds)
    else
      validationInfo <- list("split", permutations, percent)
  } else if(validation == "leaveOut") {
    validationInfo <- list("leave", leave)
  } else {
    validationInfo <- list("fold", folds)
  }

  resultsByVariety <- lapply(resultsByVariety, function(result) {if(length(unlist(result[["tune"]])) == 0) result[["tune"]] <- list(NULL); result;}) # Shorten it.
  resultsByVariety <- lapply(resultsByVariety, function(result) {if(length(unlist(result[["ranked"]])) == 0) result[["ranked"]] <- list(NULL); result;}) # Shorten it.

  classifyResults <- lapply(varietyNames, function(variety)
  {
    # Might be NULL if selection is done within training.
    selectionName <- ifelse(is.null(selectParams), "Unspecified", selectParams@selectionName)
    
    ClassifyResult(datasetName, classificationName, selectionName, rownames(measurements), allFeatures, consideredFeatures,
                   resultsByVariety[[variety]][["ranked"]], resultsByVariety[[variety]][["selected"]], resultsByVariety[[variety]][["models"]], resultsByVariety[[variety]][["predictions"]],
                   classes, validationInfo, resultsByVariety[[variety]][["tune"]])
  })

  names(classifyResults) <- varietyNames
  if(multipleVarieties == FALSE) classifyResults <- classifyResults[[1]]
  classifyResults  
})

setMethod("runTests", c("MultiAssayExperiment"),
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, restrict = NULL)
  runTests(tablesAndClasses[["dataTable"]], tablesAndClasses[["classes"]], ...)            
})

setGeneric("runTestsEasyHard", function(measurements, ...)
           {standardGeneric("runTestsEasyHard")})

setMethod("runTestsEasyHard", c("MultiAssayExperiment"),
          function(measurements, easyDatasetID = "clinical", hardDatasetID = names(measurements)[1],
                   featureSets = NULL, metaFeatures = NULL, minimumOverlapPercent = 80,
                   datasetName = NULL, classificationName = "Easy-Hard Classifier", 
                   validation = c("permute", "leaveOut", "fold"),
                   permutePartition = c("fold", "split"),
                   permutations = 100, percent = 25, folds = 5, leave = 2,
                   seed, parallelParams = bpparam(), ..., verbose = 1)
          {
            validation <- match.arg(validation)
            permutePartition <- match.arg(permutePartition)
            if(!missing(seed)) set.seed(seed)
            resultTypes <- c("selected", "models", "testSet", "predictions", "tune")
            # Elements of list returned by runTest.
            
            if(easyDatasetID == "clinical")
            {
              easyDataset <- MultiAssayExperiment::colData(measurements) # Will be DataFrame
              easyDataset <- easyDataset[!is.na(easyDataset[, "class"]), ]
              easyDataset <- easyDataset[, -match("class", colnames(easyDataset))] # Don't let the class variable go into the classifier training!
            } else if(easyDatasetID %in% names(measurements))
            {
              easyDataset <- measurements[, , easyDatasetID][[1]] # Get the underlying data container e.g. matrix.
              if(is.matrix(easyDataset))
                easyDataset <- t(easyDataset) # Make the variables be in columns.
            } else {
              stop("'easyDatasetID' is not \"clinical\" nor the name of any assay in 'measurements'.")
            }
            
            if(hardDatasetID %in% names(measurements))
            {
              hardDataset <- measurements[, , hardDatasetID][[1]] # Get the underlying data container e.g. matrix.
              hardDataset <- S4Vectors::DataFrame(t(hardDataset), check.names = FALSE) # Variables as columns.
            } else {
              stop("'hardDatasetID' is not the name of any assay in 'measurements'.")
            }
            
            # Avoid samples in one dataset but absent from the other.
            commonSamples <- intersect(rownames(easyDataset), rownames(hardDataset))
            measurements <- measurements[ , commonSamples, ]
            easyDataset <- easyDataset[commonSamples, ]
            hardDataset <- hardDataset[commonSamples, ]
            
            allFeatures <- S4Vectors::DataFrame(dataset = easyDatasetID, feature = colnames(easyDataset))
            allFeatures <- rbind(allFeatures, S4Vectors::DataFrame(dataset = hardDatasetID, feature = colnames(hardDataset)))
            classes <- MultiAssayExperiment::colData(measurements)[commonSamples, "class"]
            
            # Could refer to features or feature sets, depending on if a selection method utilising feature sets is used.
            easyFeaturesNumber <- sum(allFeatures[, "dataset"] == easyDatasetID)
            if(!is.null(featureSets))
              consideredFeatures <- length(featureSets@sets) + easyFeaturesNumber
            else
              consideredFeatures <- ncol(hardDataset) + easyFeaturesNumber
            consideredSamples <- nrow(MultiAssayExperiment::colData(measurements))
            
            if(!is.null(featureSets))
            {
              hardFeatureNames <- rownames(hardDataset)
              
              # Filter out the edges or the features from the sets which are not in measurements.
              featureSetsList <- featureSets@sets
              if(class(featureSetsList[[1]]) == "matrix")
              {
                setsSizes <- sapply(featureSetsList, nrow)
                edgesAll <- do.call(rbind, featureSetsList)
                networkNames <- rep(names(featureSetsList), setsSizes)
                edgesKeep <- edgesAll[, 1] %in% hardFeatureNames & edgesAll[, 2] %in% hardFeatureNames
                edgesFiltered <- edgesAll[edgesKeep, ]
                networkNamesFiltered <- networkNames[edgesKeep]
                setsRows <- split(1:nrow(edgesFiltered), factor(networkNamesFiltered, levels = featureSetsList))
                featureSetsListFiltered <- lapply(setsRows, function(setRows) edgesFiltered[setRows, , drop = FALSE])
                setsSizesFiltered <- sapply(featureSetsListFiltered, nrow)
              } else { # A set of features without edges, such as a gene set.
                setsSizes <- sapply(setsNodes, length)
                nodesVector <- unlist(featureSetsList)
                setsVector <- rep(names(featureSetsList), setsSizes)
                keepNodes <- !is.na(match(nodesVector, hardFeatureNames))
                nodesVector <- nodesVector[keepNodes]
                setsVector <- setsVector[keepNodes]
                featureSetsListFiltered <- split(nodesVector, factor(setsVector, levels = names(featureSetsList)))
                setsSizesFiltered <- sapply(featureSetsListFiltered, length)
              }
              keepSets <- setsSizesFiltered / setsSizes * 100 >= minimumOverlapPercent
              featureSetsListFiltered <- featureSetsListFiltered[keepSets]
              featureSets <- FeatureSetCollection(featureSetsListFiltered)
              hardDataset <- hardDataset[, hardFeatureNames %in% unlist(featureSetsListFiltered)]
              
              if(verbose >= 1 && is.null(.iteration)) # Being used by the user, not called by runTests.
                message("After filtering features, ", length(featureSetsListFiltered), " out of ", length(featureSetsList), " sets remain.")
            }
            
            if(validation == "permute")
            {
              if(permutePartition == "fold")
              {
                samplesFolds <- lapply(1:permutations, function(permutation)
                {
                  classesFolds <- lapply(levels(classes), function(className)
                  {
                    whichSamples <- which(classes == className)
                    split(sample(whichSamples), rep(1:folds, length.out = length(whichSamples)))
                  })
                  allFolds <- lapply(1:folds, function(fold)
                  {
                    unlist(lapply(classesFolds, "[[", fold))
                  })
                })
              } else { # Is split.
                samplesTrain <- (100 - percent) / 100 * table(classes)
                samplesFolds <- lapply(1:permutations, function(permutation)
                {
                  trainSet <- unlist(mapply(function(className, number)
                  {
                    sample(which(classes == className), number)
                  }, levels(classes), samplesTrain))
                  testSet <- setdiff(1:length(classes), trainSet)
                  list(trainSet, testSet)
                })
              }
              
              results <- bpmapply(function(sampleFolds, sampleNumber, ...)
              {print(sampleNumber)
                if(verbose >= 1 && sampleNumber %% 10 == 0)
                  message("Processing sample set ", sampleNumber, '.')
                if(permutePartition == "fold")
                {
                  lapply(1:folds, function(foldIndex)
                  {print(foldIndex)
                    #if(sampleNumber == 4 && foldIndex == 5) browser()
                    runTestEasyHard(measurements, easyDatasetID, hardDatasetID, featureSets, metaFeatures, minimumOverlapPercent,
                                    datasetName, classificationName,
                                    unlist(sampleFolds[-foldIndex]), sampleFolds[[foldIndex]], ..., verbose = verbose,
                                    .iteration = c(sampleNumber, foldIndex))
                  })
                } else { # Split mode.
                  runTestEasyHard(measurements, easyDatasetID, hardDatasetID, featureSets, metaFeatures, minimumOverlapPercent,
                                  datasetName, classificationName,
                                  sampleFolds[[1]], sampleFolds[[2]], ..., verbose = verbose, .iteration = sampleNumber)
                }
              }, samplesFolds, as.list(1:permutations), MoreArgs = list(...), BPPARAM = parallelParams, SIMPLIFY = FALSE)
            } else if(validation == "leaveOut") # leave k out.
            {
              testSamples <- as.data.frame(utils::combn(consideredSamples, leave))
              trainingSamples <- lapply(testSamples, function(sample) setdiff(1:consideredSamples, sample))
              
              results <- bpmapply(function(trainingSample, testSample, sampleNumber, ...)
              {
                if(verbose >= 1 && sampleNumber %% 10 == 0)
                  message("Processing sample set ", sampleNumber, '.')
                
                runTestEasyHard(measurements, easyDatasetID, hardDatasetID, featureSets, metaFeatures, minimumOverlapPercent,
                                datasetName, classificationName,
                                trainingSample, testSample, ..., verbose = verbose, .iteration = sampleNumber)
              }, trainingSamples, testSamples, (1:length(trainingSamples)), MoreArgs = list(...),
              BPPARAM = parallelParams, SIMPLIFY = FALSE)
            } else { # Unresampled, ordinary k-fold cross-validation.
              classesFolds <- lapply(levels(classes), function(className)
              {
                whichSamples <- which(classes == className)
                split(sample(whichSamples), rep(1:folds, length.out = length(whichSamples)))
              })
              samplesFolds <- lapply(1:folds, function(fold)
              {
                unlist(lapply(classesFolds, "[[", fold))
              })
              
              if(verbose >= 1)
                message("Processing ", folds, "-fold cross-validation.")
              results <- bplapply(1:folds, function(foldIndex, ...)
              {
                runTestEasyHard(measurements, easyDatasetID, hardDatasetID, featureSets, metaFeatures, minimumOverlapPercent,
                                datasetName, classificationName,
                                unlist(samplesFolds[-foldIndex]), samplesFolds[[foldIndex]], ..., verbose = verbose,
                        .iteration = foldIndex)
              }, measurements, ..., BPPARAM = parallelParams)
            }
            
            if(validation == "permute" && permutePartition == "split" || validation %in% c("leaveOut", "fold"))
            {
              resultErrors <- sapply(results, function(result) is.character(result))
              if(sum(resultErrors) == length(results))
              {
                message("Error: All cross-validations had an error.")
	        if(length(unique(unlist(results))) == 1)
                  message("The common problem is: ", unlist(results)[[1]])
                return(results)
              } else if(sum(resultErrors) != 0) # Filter out cross-validations resulting in error.
              {
		warning(paste(sum(resultErrors),  "cross-validations had an error and have been removed from the results."))
                results <- results[!resultErrors]
              }
            } else { # Result has nested lists, because of permutation reordering with folding.
              resultErrors <- lapply(results, function(resample) lapply(resample, is.character))
              if(sum(unlist(resultErrors)) == permutations * folds)
              {
                message("Error: All cross-validations had an error.")
	      	if(length(unique(unlist(results))) == 1)
                  message("The common problem is: ", unlist(results)[[1]])
                return(results)
              } else if(sum(unlist(resultErrors)) != 0) # Filter out error cross-validations.
              {
		warning(paste(sum(unlist(resultErrors)),  "cross-validations had an error and have been removed from the results."))
                results <- results[sapply(results, function(resample) !all(sapply(resample, is.character)))]
                results <- lapply(results, function(resample) resample[!sapply(resample, is.character)])
              }
            }
            
            if(validation == "permute")
            {
              if(permutePartition == "fold")
              {
                if(class(results[[1]][[1]][["predictions"]]) == "list") 
                {
                  multipleVarieties <- TRUE
                  varietyNames <- names(results[[1]][[1]][["predictions"]])
                } else {
                  multipleVarieties <- FALSE
                  varietyNames <- "None"
                }
              } else {
                if(class(results[[1]][["predictions"]]) == "list")
                {
                  multipleVarieties <- TRUE
                  varietyNames <- names(results[[1]][["predictions"]])
                } else {
                  multipleVarieties <- FALSE
                  varietyNames <- "None"
                }
              }
              
              resultsByVariety <- lapply(varietyNames, function(varietyName)
              {
                resultsLists <- lapply(resultTypes, function(resultType)
                {
                  resultList <- lapply(results, function(sample)
                  {
                    if(resultType == "testSet" || multipleVarieties == FALSE)
                    {
                      if(resultType %in% c("testSet", "predictions"))
                      {
                        if(permutePartition == "fold")
                        {
                          if(resultType == "testSet" || class(sample[[1]][["predictions"]]) != "data.frame")
                          {
                            # First fold predictions are a not data.frame.
                            reshaped <- unlist(lapply(sample, function(fold) fold[[resultType]]))
                            if(resultType == "predictions")
                              reshaped <- factor(reshaped, levels = levels(classes))
                            reshaped
                          } else { # Predictions is a data.frame with classes and scores.
                            resultTable <- do.call(rbind, lapply(sample, function(fold) fold[[resultType]]))
                            resultTable[, "class"] <- factor(resultTable[, "class"], levels = levels(classes))
                            resultTable
                          }
                        } else { # permutePartition is "split".
                          information <- sample[[resultType]]
                          if(resultType == "predictions")
                            information <- factor(information, levels = levels(classes))
                          information
                        }
                      } else {
                        if(permutePartition == "fold")
                          lapply(sample, function(fold) fold[[resultType]])
                        else
                          sample[[resultType]]
                      }
                    } else { # There are multiple varieties.
                      if(resultType == "predictions") # Predictions with varieties.
                      {
                        if(permutePartition == "fold")
                        {
                          prediction <- lapply(sample, function(fold) fold[[resultType]][[varietyName]])
                          if(class(prediction[[1]]) == "data.frame")
                          {
                            reshaped <- do.call(rbind, prediction)
                            reshaped[, "class"] <- factor(reshaped[, "class"], levels = levels(classes))
                            reshaped
                          } else
                            factor(unlist(prediction), levels = levels(classes))
                        } else { # permutePartition is "split".
                          information <- sample[[resultType]][[varietyName]]
                          if(resultType == "predictions")
                            information <- factor(information, levels = levels(classes))
                          information                
                        }
                      } else { # Ranked, selected features, models or test samples, with varieties.
                        if(permutePartition == "fold")
                        {
                          if(is.null(sample[[1]][[resultType]]))
                            NULL
                          else
                            lapply(sample, function(fold) fold[[resultType]][[varietyName]])
                        } else {
                          if(is.null(sample[[resultType]]))
                            NULL
                          else
                            sample[[resultType]][[varietyName]]
                        }
                      }
                    }
                  })
                  if(is.null(unlist(resultList, use.names = FALSE))) list(NULL) else resultList # Collapse multiple NULLs into one NULL.
                })
                names(resultsLists) <- resultTypes
                resultsLists
              })
              names(resultsByVariety) <- varietyNames
            } else {# leave k out or ordinary, unresampled k-fold cross-validation.
              if(class(results[[1]][["predictions"]]) == "list")
              {
                multipleVarieties <- TRUE
                varietyNames <- names(results[[1]][["predictions"]])
              } else {
                multipleVarieties <- FALSE
                varietyNames <- "None"
              }
              
              resultsByVariety <- lapply(varietyNames, function(varietyName)
              {
                resultsLists <- lapply(resultTypes, function(resultType)
                {
                  resultList <- lapply(results, function(sample) sample[[resultType]])
                  if(multipleVarieties == TRUE && resultType != "testSet") # Special majority class classification. Expand to number of varieties so that loops below work.
                  {
                    whichNotVariety <- which(sapply(resultList, function(result) !class(result) == "list"))
                    resultList[whichNotVariety] <- lapply(resultList[whichNotVariety], function(result) setNames(rep(list(result), length(varietyNames)), varietyNames))
                  }
                  if(resultType %in% c("testSet", "predictions") && multipleVarieties == FALSE && class(resultList[[1]]) != "data.frame")
                  {
                    reshaped <- unlist(resultList)
                    if(resultType == "predictions")
                      reshaped <- factor(reshaped, levels = levels(classes)) # Enforce that levels are in the same order as input by the user.
                    reshaped
                  } else if (resultType == "predictions" && multipleVarieties == FALSE) { # Predictions are in a data.frame.
                    resultTable <- do.call(rbind, resultList)
                    resultTable[, "class"] <- factor(resultTable[, "class"], levels = levels(classes))
                    resultTable
                  } else if(multipleVarieties == TRUE && resultType == "predictions" && class(resultList[[1]]) != "data.frame") {
                    factor(unlist(lapply(resultList, "[[", varietyName)), levels = levels(classes))
                  } else if(multipleVarieties == TRUE && resultType == "predictions") { # Predictions are factor or numeric.
                    resultTable <- do.call(rbind, lapply(resultList, "[[", varietyName))
                    resultTable[, "class"] <- factor(resultTable[, "class"], levels = levels(classes))
                    resultTable
                  } else if(multipleVarieties == TRUE && resultType %in% c("ranked", "selected", "tune", "models")) {
                    returnValue <- lapply(resultList, "[[", varietyName)
                    if(is.null(unlist(returnValue, use.names = FALSE))) NULL else returnValue
                  } else { # For multipleVarieties being FALSE, and not the predicted classes.
                    resultList
                  }
                })
                names(resultsLists) <- resultTypes
                resultsLists
              })
              names(resultsByVariety) <- varietyNames
            }
            
            predictionTablesByVariety <- lapply(resultsByVariety, function(resultVariety)
            {
              if(validation == "permute")
              {
                lapply(1:length(results), function(resample)
                {
                  sampleNames <- unlist(resultVariety[["testSet"]][[resample]])
                  switch(class(resultVariety[["predictions"]][[resample]]),
                         factor = data.frame(sample = sampleNames, class = factor(resultVariety[["predictions"]][[resample]], levels = levels(classes)), stringsAsFactors = FALSE, row.names = NULL),
                         numeric = data.frame(sample = sampleNames, resultVariety[["predictions"]][[resample]], stringsAsFactors = FALSE, row.names = NULL),
                         data.frame = data.frame(sample = sampleNames,
                                                 class = factor(resultVariety[["predictions"]][[resample]][, sapply(resultVariety[["predictions"]][[resample]], class) == "factor"], levels = levels(classes)),
                                                 resultVariety[["predictions"]][[resample]][, sapply(resultVariety[["predictions"]][[resample]], class) == "numeric"], stringsAsFactors = FALSE, row.names = NULL))
                  
                })
              } else { # leave k out or ordinary, unresampled k-fold cross-validation.
                sampleNames <- unlist(resultVariety[["testSet"]])
                list(switch(class(resultVariety[["predictions"]]),
                            factor = data.frame(sample = sampleNames, class = factor(resultVariety[["predictions"]]), stringsAsFactors = FALSE, row.names = NULL),
                            numeric = data.frame(sample = sampleNames, resultVariety[["predictions"]], stringsAsFactors = FALSE, row.names = NULL),
                            data.frame = data.frame(sample = sampleNames,
                                                    class = factor(resultVariety[["predictions"]][, sapply(resultVariety[["predictions"]], class) == "factor"], levels = levels(classes)),
                                                    resultVariety[["predictions"]][, sapply(resultVariety[["predictions"]], class) == "numeric"], stringsAsFactors = FALSE, row.names = NULL)))
              }
            })
            
            resultsByVariety <- mapply(function(resultVariety, predictionVariety) # Put formatted prediction tables into results list.
            {
              resultVariety[["predictions"]] <- predictionVariety
              resultVariety
            }, resultsByVariety, predictionTablesByVariety, SIMPLIFY = FALSE)
            
            if(validation == "permute")
            {
              if(permutePartition == "fold")
                validationInfo <- list("permuteFold", permutations, folds)
              else
                validationInfo <- list("split", permutations, percent)
            } else if(validation == "leaveOut") {
              validationInfo <- list("leave", leave)
            } else {
              validationInfo <- list("fold", folds)
            }
            
            resultsByVariety <- lapply(resultsByVariety, function(result) {if(length(unlist(result[["tune"]])) == 0) result[["tune"]] <- list(NULL); result;}) # Shorten it.
            resultsByVariety <- lapply(resultsByVariety, function(result) {if(length(unlist(result[["ranked"]])) == 0) result[["ranked"]] <- list(NULL); result;}) # Shorten it.
            
            classifyResults <- lapply(varietyNames, function(variety)
            {
              selectionName <- "Sample Grouping Purity for Easy Data Set"
              selectParams <- list(...)[["hardClassifierParams"]]
              whichSelect <- which(sapply(selectParams, class) == "SelectParams")
              
              if(length(whichSelect) > 0)
              {
                selectParams <- selectParams[[whichSelect]]
                selectionName <- paste(selectionName, "and", selectParams@selectionName, "for Hard Data Set")
              } else {
                selectionName <- paste(selectionName, '.', sep = '')
              }
              
              ClassifyResult(datasetName, classificationName, selectionName, rownames(MultiAssayExperiment::colData(measurements)), allFeatures, consideredFeatures,
                             resultsByVariety[[variety]][["ranked"]], resultsByVariety[[variety]][["selected"]], resultsByVariety[[variety]][["models"]], resultsByVariety[[variety]][["predictions"]],
                             classes, validationInfo, resultsByVariety[[variety]][["tune"]])
            })
            
            names(classifyResults) <- varietyNames
            if(multipleVarieties == FALSE) classifyResults <- classifyResults[[1]]
            classifyResults  
          })

