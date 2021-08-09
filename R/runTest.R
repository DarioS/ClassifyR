setGeneric("runTest", function(measurements, ...)
           standardGeneric("runTest"))

setMethod("runTest", c("matrix"), # Matrix of numeric measurements.
  function(measurements, classes, ...)
{
  if(is.null(colnames(measurements)))
    stop("'measurements' matrix must have sample identifiers as its column names.")    
  runTest(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("runTest", c("DataFrame"), # Clinical data or one of the other inputs, transformed..
function(measurements, classes,
         featureSets = NULL, metaFeatures = NULL, minimumOverlapPercent = 80,
         datasetName, classificationName, training, testing,
         params = list(SelectParams(), TrainParams(), PredictParams()),
         verbose = 1, .iteration = NULL)
{
  if(is.null(rownames(measurements)))
    stop("'measurements' DataFrame must have sample identifiers as its row names.")  
  splitDataset <- .splitDataAndClasses(measurements, classes)

  stagesParamClasses <- sapply(params, class)
  if(match("TrainParams", stagesParamClasses) > match("PredictParams", stagesParamClasses))
    stop("\"PredictParams\" variable must not be before \"TrainParams\" in 'params' list.")

  transformParams <- params[[match("TransformParams", stagesParamClasses)]]
  selectParams <- params[[match("SelectParams", stagesParamClasses)]]
  trainParams <- params[[match("TrainParams", stagesParamClasses)]]
  predictParams <- params[[match("PredictParams", stagesParamClasses)]]
  
  # All input features.
  if(!is.null(S4Vectors::mcols(measurements)))
    allFeatures <- S4Vectors::mcols(measurements)
  else
    allFeatures <- colnames(measurements)
  
  # Could refer to features or feature sets, depending on if a selection method utilising feature sets is used.
  if(!is.null(selectParams) && !is.null(featureSets))
    consideredFeatures <- length(featureSets@sets)
  else
    consideredFeatures <- ncol(measurements)
  
  if(!is.null(featureSets) && is.null(.iteration)) # Feature sets provided and runTest is being called by the user, so need to be filtered now.
  {
    if(!is.null(S4Vectors::mcols(measurements)))
      featureNames <- S4Vectors::mcols(measurements)[, "feature"]
    else
      featureNames <- colnames(measurements)
    
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
      setsRows <- split(1:nrow(edgesFiltered), factor(networkNamesFiltered, levels = featureSetsList))
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
    
    if(verbose >= 1 && is.null(.iteration)) # Being used by the user, not called by runTests.
      message("After filtering features, ", length(featureSetsListFiltered), " out of ", length(featureSetsList), " sets remain.")
  }
  
  lastSize <- 1
  for(stageIndex in 1:length(params))
  {
    switch(stagesParamClasses[[stageIndex]],
           TransformParams = {
                               if(length(transformParams@intermediate) != 0)
                               {
                                 intermediates <- mget(transformParams@intermediate)
                                  if(!is.null(names(transformParams@intermediate)))
                                    names(intermediates) <- names(transformParams@intermediate)
                                  transformParams@otherParams <- c(transformParams@otherParams, intermediates)                                 
                               }

                               measurements <- tryCatch(.doTransform(measurements, training, transformParams, verbose), error = function(error) error[["message"]])
                               if(is.character(measurements)) return(measurements) # An error occurred.
                               newSize <- if(class(measurements) == "list") length(measurements) else 1
                             },
              SelectParams = {
                               if(length(selectParams@intermediate) != 0)
                               {
                                 intermediates <- mget(selectParams@intermediate)
                                  if(!is.null(names(selectParams@intermediate)))
                                    names(intermediates) <- names(selectParams@intermediate)
                                  selectParams@otherParams <- c(selectParams@otherParams, intermediates)
                               }

                               topFeatures <- tryCatch(.doSelection(measurements, classes, featureSets, metaFeatures, training, selectParams,
                                                                trainParams, predictParams, verbose), error = function(error) error[["message"]])
                               if(is.character(topFeatures)) return(topFeatures) # An error occurred.
  
                               if(class(topFeatures[[2]]) == "list") # Check the chosen features list element, because a ranking is not present for ensemble selection.
                               {
                                 multiSelection <- TRUE
                               } else {
                                 multiSelection <- FALSE
                               }

                               rankedFeatures <- topFeatures[[1]] # Extract for result object.
                               selectedFeatures <- topFeatures[[2]] # Extract for subsetting.

                               if(selectParams@subsetToSelections == TRUE)
                               {
                                 if(multiSelection == FALSE)
                                 {
                                   if(is.null(metaFeatures))
                                   {
                                     if(class(measurements) != "list") # Put into list.
                                       measurements <- list(measurements)
                                     measurements <- lapply(measurements, function(variety)
                                                     {
                                                       if(is.null(S4Vectors::mcols(variety)) == TRUE)
                                                       { # Input was ordinary matrix or DataFrame and no network features were used.
                                                         variety[, selectedFeatures, drop = FALSE]
                                                       } else { # Input was MultiAssayExperiment.
                                                         selectedColumns <- apply(selectedFeatures, 1, function(selectedFeature)
                                                         {
                                                           intersect(which(selectedFeature[1] == S4Vectors::mcols(variety)[, "dataset"]),
                                                                     which(selectedFeature[2] == S4Vectors::mcols(variety)[, "feature"]))
                                                         })
                                                         variety <- variety[, selectedColumns, drop = FALSE]
                                                         variety
                                                       }
                                                     })
                                     if(length(measurements) == 1 && class(measurements) == "list")  # Restore to original container type.
                                       measurements <- measurements[[1]]
                                     measurements
                                   } else {
                                     metaFeatures <- metaFeatures[, S4Vectors::mcols(metaFeatures)[, "original"] %in% selectedFeatures]
                                   }
                                 } else { # Multiple varieties of selections.
                                   if(is.null(metaFeatures))
                                   {
                                     if(class(measurements) != "list") # Put into list.
                                       measurements <- list(measurements)
                                     
                                     measurements <- lapply(measurements, function(variety)
                                                     {
                                                       lapply(selectedFeatures, function(features)
                                                       {
                                                           if(is.null(S4Vectors::mcols(variety)) == TRUE)
                                                           { # Input was ordinary matrix or DataFrame.
                                                             variety[, features, drop = FALSE]
                                                           } else { # Input was MultiAssayExperiment.
                                                             selectedColumns <- apply(selectedFeatures, 2, function(selectedFeature)
                                                             {
                                                               intersect(which(selectedFeature[1] == S4Vectors::mcols(variety)[, "dataset"]),
                                                                         which(selectedFeature[2] == S4Vectors::mcols(variety)[, "feature"]))
                                                             })
                                                             variety <- variety[, selectedColumns, drop = FALSE]
                                                             variety
                                                           }
                                                         })
                                                       })
                                     if(length(measurements) == 1 && class(measurements) == "list")  # Restore to original container type.
                                       measurements <- measurements[[1]]
                                   } else {
                                     metaFeatures <- lapply(selectedFeatures, function(features)
                                                     {
                                                       metaFeatures[, S4Vectors::mcols(metaFeatures)[, "original"] %in% features]
                                                     })
                                   }
                                 }
                               } else { # Don't subset to the selected features.
                                 if(multiSelection == TRUE) # Multiple selection varieties. Replicate the experimental data.
                                 {
                                   if(is.null(metaFeatures))
                                   {
                                     if(class(measurements) != "list")
                                       measurements <- lapply(selectedFeatures, function(features) measurements)
                                     else
                                       measurements <- lapply(measurements, function(variety)
                                                            lapply(selectedFeatures, function(features) variety))
                                   } else {
                                     if(class(metaFeatures) != "list")
                                       metaFeatures <- lapply(selectedFeatures, function(features) metaFeatures)
                                     else
                                       metaFeatures <- lapply(metaFeatures, function(variety)
                                                            lapply(selectedFeatures, function(features) variety))
                                   }
                                 }
                               }

                               if(is.null(metaFeatures))
                               {
                                 if(class(measurements) == "list" && class(measurements[[1]]) == "list")
                                 {
                                   oldNames <- sapply(measurements, names)
                                   newNames <- unlist(lapply(measurements, names))
                                   measurements <- unlist(measurements, recursive = FALSE)
                                   names(measurements) <- paste(rep(oldNames, each = length(measurements[[1]])), newNames, sep = ',')
                                 }
                                 
                                 if(class(measurements) == "list") newSize <- length(measurements) else newSize <- 1
                                 lastSize <- newSize
                               } else {
                                 if(class(metaFeatures) == "list" && class(metaFeatures[[1]]) == "list")
                                 {
                                   oldNames <- sapply(metaFeatures, names)
                                   newNames <- unlist(lapply(metaFeatures, names))
                                   metaFeatures <- unlist(metaFeatures, recursive = FALSE)
                                   names(metaFeatures) <- paste(rep(oldNames, each = length(metaFeatures[[1]])), newNames, sep = ',')
                                 }
                                 
                                 if(class(metaFeatures) == "list") newSize <- length(metaFeatures) else newSize <- 1
                                 lastSize <- newSize                                 
                               }
                             }, 
              TrainParams = {
                              if(length(trainParams@intermediate) != 0)
                              {
                                intermediates <- mget(trainParams@intermediate)
                                if(!is.null(names(trainParams@intermediate)))
                                  names(intermediates) <- names(trainParams@intermediate)
                                trainParams@otherParams <- c(trainParams@otherParams, intermediates)
                              }

                              if(is.null(metaFeatures))
                                useData <- measurements
                              else # Used some derived features instead.
                                useData <- metaFeatures
                              trained <- tryCatch(.doTrain(useData, classes, training, testing, trainParams, predictParams, verbose),
                                                  error = function(error) error[["message"]])
                              if(is.character(trained)) return(trained) # An error occurred.

                              newSize <- if("list" %in% class(trained)) length(trained) else 1
                              if(newSize / lastSize != 1) # More varieties were created.
                              {
                                if(is.null(metaFeatures))
                                {
                                measurements <- unlist(lapply(if(class(measurements) == "list") measurements else list(measurements), function(variety)
                                                                                      lapply(1:(newSize / lastSize), function(replicate) variety)),
                                                                               recursive = FALSE)
                                names(measurements) <- names(trained)
                                } else {
                                  metaFeatures <- unlist(lapply(if(class(metaFeatures) == "list") metaFeatures else list(metaFeatures), function(variety)
                                                                                      lapply(1:(newSize / lastSize), function(replicate) variety)),
                                                                               recursive = FALSE)
                                  names(metaFeatures) <- names(trained)
                                }
                              }

                              lastSize <- newSize
                              if("list" %in% class(trained))
                              {
                                tuneDetails <- lapply(trained, attr, "tune")
                                if(!is.null(trainParams@getFeatures)) # Features chosen inside classifier.
                                {
                                  featureInfo <- lapply(trained, trainParams@getFeatures)
                                  rankedFeatures <- lapply(featureInfo, '[[', 1)
                                  selectedFeatures <- lapply(featureInfo, '[[', 2)
                                }
                              } else {
                                tuneDetails <- attr(trained, "tune")
                                if(!is.null(trainParams@getFeatures)) # Features chosen inside classifier.
                                {                                
                                  rankedFeatures <- trainParams@getFeatures(trained)[[1]]
                                  selectedFeatures <- trainParams@getFeatures(trained)[[2]]
                                }
                              }
                                if(is.null(tuneDetails)) tuneDetails <- list(tuneDetails)
                              },
              PredictParams = {
                                if(is.null(metaFeatures))
                                  useData <- measurements
                                else # Used some derived features instead.
                                  useData <- metaFeatures
                                
                                if(length(predictParams@intermediate) != 0)
                                {
                                  intermediates <- mget(predictParams@intermediate)
                                  if(!is.null(names(predictParams@intermediate)))
                                    names(intermediates) <- names(predictParams@intermediate)
                                  predictParams@otherParams <- c(predictParams@otherParams, intermediates)                                  
                                }
                                 
                                predictedClasses <- tryCatch(.doTest(trained, useData, testing, predictParams, verbose),
                                                             error = function(error) error[["message"]])
                                if(is.character(predictedClasses)) # An error occurred.
                                  return(predictedClasses) # Return early. Don't make a ClassifyResult below.
                              }
           )
    
  }
  if(is.logical(testing)) testing <- which(testing)
  if(is.numeric(testing)) testing <- rownames(measurements)[testing]
  # Rankings and selections might not be explicitly returned, such as for random forest classifier.
  if(!exists("rankedFeatures")) rankedFeatures <- NULL
  if(!exists("selectedFeatures")) selectedFeatures <- NULL
  if(is.null(predictParams@predictor)) models <- NULL else models <- trained # One function for training and testing. Typically, the models aren't returned to the user, such as Poisson LDA implemented by PoiClaClu.
  if(!is.null(.iteration)) # This function was called by runTests.
  {
    list(ranked = rankedFeatures, selected = selectedFeatures, models = models, testSet = testing, predictions = predictedClasses, tune = tuneDetails)
  } else { # runTest is being used directly, rather than from runTests. Create a ClassifyResult object.
    if(class(predictedClasses) != "list")
    {
      return(ClassifyResult(datasetName, classificationName, selectParams@selectionName, rownames(measurements), allFeatures, consideredFeatures,
                            list(rankedFeatures), list(selectedFeatures), list(models), list(data.frame(sample = testing, class = predictedClasses)),
                            classes, list("independent"), tuneDetails)
             )
    } else { # A variety of predictions were made.
      if(!"list" %in% class(selectedFeatures))
      {
        rankedFeatures <- list(rankedFeatures)
        selectedFeatures <- list(selectedFeatures)
        models <- list(models)
      }
      return(mapply(function(varietyPredictions, varietyTunes, varietyRanked, varietySelected, varietyModels)
      {
        if(is.null(varietyTunes)) varietyTunes <- list(varietyTunes)
        ClassifyResult(datasetName, classificationName, selectParams@selectionName, rownames(measurements), allFeatures, consideredFeatures,
                       list(varietyRanked), list(varietySelected), list(varietyModels), list(data.frame(sample = testing, class = varietyPredictions)),
                       classes, list("independent"), varietyTunes)
      }, predictedClasses, tuneDetails, rankedFeatures, selectedFeatures, models, SIMPLIFY = FALSE))
    }
  }  
})

setMethod("runTest", c("MultiAssayExperiment"),
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, restrict = NULL)
  runTest(tablesAndClasses[["dataTable"]], tablesAndClasses[["classes"]], ...)            
})

setGeneric("runTestEasyHard", function(measurements, ...)
standardGeneric("runTestEasyHard"))

setMethod("runTestEasyHard", c("MultiAssayExperiment"),
          function(measurements, easyDatasetID = "clinical", hardDatasetID = names(measurements)[1],
                   featureSets = NULL, metaFeatures = NULL, minimumOverlapPercent = 80,
                   datasetName = NULL, classificationName = "Easy-Hard Classifier", training, testing, ..., verbose = 1, .iteration = NULL)
          {
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
            classes <- MultiAssayExperiment::colData(measurements)[, "class"]
            
            # Could refer to features or feature sets, depending on if a selection method utilising feature sets is used.
            easyFeaturesNumber <- sum(allFeatures[, "dataset"] == easyDatasetID)
            if(!is.null(featureSets))
              consideredFeatures <- length(featureSets@sets) + easyFeaturesNumber
            else
              consideredFeatures <- ncol(hardDataset) + easyFeaturesNumber
            
            if(!is.null(featureSets) && is.null(.iteration)) # Feature sets provided and runTestEasyHard is being called by the user, so need to be filtered now.
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

            trained <- easyHardClassifierTrain(measurements[ , training, ], easyDatasetID, hardDatasetID, featureSets, metaFeatures, minimumOverlapPercent, NULL, classificationName, ..., verbose = verbose)
            hardParams <- list(...)[["hardClassifierParams"]]
            if(is.null(hardParams)) # They were not specified by the user. Use default value.
              predictParams <- PredictParams()
            else
              predictParams <- hardParams[[which(sapply(hardParams, class) == "PredictParams")]]
            test <- measurements[ , testing, ]
            trainClass <- class(trained)
            if(trainClass == "EasyHardClassifier")
            {
              trained <- list(trained)
              predictedClasses <- lapply(trained, function(model) easyHardClassifierPredict(model, test, predictParams, verbose))
            } else { # Modify the predictParams' otherParams.
              trainedVarietyParams <- strsplit(names(trained), "=|,")
              predictParams <- lapply(trainedVarietyParams, function(varietyParams)
                               {
                                 combinationParams <- predictParams
                                 for(paramIndex in seq(1, length(varietyParams), 2))
                                 {
                                   if(varietyParams[paramIndex] %in% names(combinationParams@otherParams))
                                   {
                                     combinationParams@otherParams[[varietyParams[paramIndex]]] <- varietyParams[paramIndex + 1]
                                   }
                                 }
                                 combinationParams
                               })
              predictedClasses <- mapply(function(varietyModel, varietyParams) easyHardClassifierPredict(varietyModel, test, varietyParams, verbose), trained, predictParams, SIMPLIFY = FALSE)
            }

            tuneDetails <- lapply(trained, attr, "tune")
            if(trainClass == "EasyHardClassifier")
            {
              trained <- trained[[1]]
              selectedFeatures <- easyHardFeatures(trained)[[2]]
              predictedClasses <- unlist(predictedClasses, recursive = FALSE)
              tuneDetails <- unlist(tuneDetails, recursive = FALSE)
            } else {
              selectedFeatures <- lapply(lapply(trained, easyHardFeatures), "[[", 2)
            }
            if(is.null(tuneDetails)) tuneDetails <- list(tuneDetails)

            if(is.logical(testing)) testing <- which(testing)
            if(is.numeric(testing)) testing <- rownames(MultiAssayExperiment::colData(measurements))[testing]
            if(!is.null(.iteration)) # This function was called by runTestsEasyHard.
            {
              list(selected = selectedFeatures, models = trained, testSet = testing, predictions = predictedClasses, tuneDetails = tuneDetails)
            } else { # runTestEasyHard is being used directly, rather than from runTestsEasyHard. Create a ClassifyResult object.
              selectionName <- "Sample Grouping Purity for Easy Data Set"
              selectParams <- list(...)[["hardClassifierParams"]]
              whichSelect <- which(sapply(selectParams, class) == "SelectParams")
              
              if(length(whichSelect) > 0)
              {
                selectParams <- selectParams[[whichSelect]]
                selectionName <- paste(selectionName, "and", selectParams@selectionName, "for Hard Data Set")
              }
              
              if(class(predictedClasses) != "list")
              {
                return(ClassifyResult(datasetName, "Easy-Hard Classifier", selectionName, rownames(MultiAssayExperiment::colData(measurements)), allFeatures, consideredFeatures,
                                      list(NULL), selectedFeatures, list(trained), list(data.frame(sample = testing, class = predictedClasses)),
                                      classes, list("independent"), tuneDetails)
                )
              } else { # A variety of predictions were made.
                if(!"list" %in% class(selectedFeatures))
                {
                  selectedFeatures <- list(selectedFeatures)
                }
                return(mapply(function(varietyPredictions, varietyTunes, varietySelected, varietyModel)
                {
                  if(is.null(varietyTunes)) varietyTunes <- list(varietyTunes)
                  ClassifyResult(datasetName, "Easy-Hard Classifier", selectionName, rownames(MultiAssayExperiment::colData(measurements)), allFeatures, consideredFeatures,
                                 list(NULL), list(varietySelected), list(varietyModel), list(data.frame(sample = testing, class = varietyPredictions)),
                                 classes, list("independent"), varietyTunes)
                }, predictedClasses, tuneDetails, selectedFeatures, trained, SIMPLIFY = FALSE))
              }
            }
          })
          