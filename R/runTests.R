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
          runTest(measurements, classes, featureSets, metaFeatures, training = unlist(sampleFolds[-foldIndex]),
                  testing = sampleFolds[[foldIndex]], params = params, verbose = verbose,
                  .iteration = c(sampleNumber, foldIndex))
        })
      } else { # Split mode.
        runTest(measurements, classes, featureSets, metaFeatures, training = sampleFolds[[1]],
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
      runTest(measurements, classes, featureSets, metaFeatures, training = trainingSample, testing = testSample,
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
      runTest(measurements, classes, featureSets, metaFeatures, training = unlist(samplesFolds[-foldIndex]),
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
      return(results)
    } else if(sum(resultErrors) != 0) # Filter out cross-validations resulting in error.
    {
      results <- results[!resultErrors]
    }
  } else { # Result has nested lists, because of permutation reordering with folding.
    resultErrors <- lapply(results, function(resample) lapply(resample, is.character))
    if(sum(unlist(resultErrors)) == permutations * folds)
    {
      message("Error: All cross-validations had an error.")
      return(results)
    } else if(sum(unlist(resultErrors)) != 0) # Filter out error cross-validations.
    {
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
        if(is.null(unlist(resultList))) list(NULL) else resultList # Collapse multiple NULLs into one NULL.
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
          if(is.null(unlist(returnValue))) NULL else returnValue
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
        sampleNames <- rownames(measurements)[resultVariety[["testSet"]][[resample]]]
        switch(class(resultVariety[["predictions"]][[resample]]),
               factor = data.frame(sample = sampleNames, class = factor(resultVariety[["predictions"]][[resample]], levels = levels(classes)), stringsAsFactors = FALSE, row.names = NULL),
               numeric = data.frame(sample = sampleNames, resultVariety[["predictions"]][[resample]], stringsAsFactors = FALSE, row.names = NULL),
               data.frame = data.frame(sample = sampleNames,
                                       class = factor(resultVariety[["predictions"]][[resample]][, sapply(resultVariety[["predictions"]][[resample]], class) == "factor"], levels = levels(classes)),
                                       resultVariety[["predictions"]][[resample]][, sapply(resultVariety[["predictions"]][[resample]], class) == "numeric"], stringsAsFactors = FALSE, row.names = NULL))
                
      })
    } else { # leave k out or ordinary, unresampled k-fold cross-validation.
      sampleNames <- rownames(measurements)[resultVariety[["testSet"]]]
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

  if(!"tuneParams" %in% names(params[[match("TrainParams", stagesParamClasses)]]@otherParams)) # No tuning specified.
    resultsByVariety <- lapply(resultsByVariety, function(result) {result[["tune"]] <- list(NULL); result;}) # Shorten it.
  resultsByVariety <- lapply(resultsByVariety, function(result) # Check for empty rankings, and shorten the list.
  {
    if(length(unlist(result[["ranked"]])) == 0)
      result[["ranked"]] <- list(NULL)
    result
  })

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
