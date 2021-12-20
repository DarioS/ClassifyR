setGeneric("quickClassify", function(measurements, ...)
  standardGeneric("quickClassify"))

setMethod("quickClassify", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
          {
            quickClassify(DataFrame(t(measurements), check.names = FALSE), classes, ...)
          })

setMethod("quickClassify", "DataFrameOrDataFrameList", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, permutations = 100, folds = 5, nFeatures = seq(10, 100, 10), cores = 1,
                   classifier = c("DLDA", "FishersLDA", "NSC", "naiveBayes", "elasticNet", "SVM", "randomForest", "kTSP"),
                   datasetMode = c("each", "combine", "allPairs"), ...)
          {
            if(is.null(rownames(measurements)))
              stop("'measurements' DataFrame must have sample identifiers as its row names.")            
            if(any(is.na(measurements)))
              stop("Some data elements are missing and classifiers don't work with missing data. Consider imputation or filtering.")
            datasetMode <- match.arg(datasetMode)
            
            if(!is.list(nFeatures)) nFeatures <- list(nFeatures) # Package it up for looping.
            if(!is.list(classes)) classes <- list(classes) # Package it up for looping.
            
            if(is(measurements, "DataFrame")) # Could have multiple data sets, details in mcols, if originating from MultiAssayExperiment.
            {
              if(!is.null(mcols(measurements)) && datasetMode == "each") datasetIDs <- unique(mcols(measurements[, "dataset"])) else datasetIDs <- "dataset"
            } else { # measurements is a DataFrameList of different projects, same assay.
              if(datasetMode != "combine") datasetIDs <- names(measurements) else datasetIDs <- "dataset"
            }
            
            if(datasetMode == "each") # Classify each data set individually.
            {
              if(length(nFeatures) < length(measurements)) nFeatures <- rep(nFeatures, length.out = length(datasetIDs))
              if(length(classifier) < length(measurements)) classifier <- rep(classifier, length.out = length(datasetIDs))
            } else if(datasetMode == "allPairs") { # All pairs of data sets.
              diffCombinations <- expand.grid(trainDatasetID = datasetIDs, predictDatasetID = datasetIDs, classifier = classifier)
              diffCombinations <- diffCombinations[diffCombinations[, "trainDatasetID"] != diffCombinations[, "predictDatasetID"], ]
            } else # Is combine.
              if(is(measurements, "DataFrameList")) { # Flatten to table of common features.
              featuresCounts <- table(unlist(lapply(measurements, colnames)))
              allCommonFeatures <- names(featuresCounts)[featuresCounts == length(featuresCounts)]
              measurements <- lapply(measurements, function(measurementsTable) measurementsTable[, allCommonFeatures, drop = FALSE])
              measurements <- list(dataset = do.call(rbind, measurements))
              classes <- list(dataset = unlist(classes))
              datasetIDs <- "dataset"
              if(length(classifier) > 1)
              {
                measurements <- rep(measurements, length.out = length(classifier))
                classes <- rep(classes, length.out = length(classifier))
                nFeatures <- rep(nFeatures, length.out = length(classifier))
                datasetIDs <- rep(datasetIDs, length.out = length(classifier))
              }
            }
            
            if(cores == 1)
            {
              BPparam <- SerialParam()
            } else { # Parallel processing is desired.
              # Also set the BPparam RNGseed if the user ran set.seed(someNumber) themselves.
              if(missing(.Random.seed)) seed <- NULL else seed <- .Random.seed[1]
              if(Sys.info()["sysname"] == "Windows") {# Only SnowParam suits Windows.
                BPparam <- SnowParam(cores, RNGseed = seed)
              } else if (Sys.info()["sysname"] %in% c("MacOS", "Linux")) {
                BPparam <- MulticoreParam(cores, RNGseed = seed) # Multicore is faster than SNOW, but it doesn't work on Windows.
              } else { # Something weird.
                BPparam <- bpparam() # BiocParallel will figure it out.
            }
            crossValParams <- CrossValParams(permutations = permutations, folds = folds, parallelParams = BPparam)
            
            data(HuRI) # Loads a variable called interactors from data folder, human reference interactome.
            # Loop for single data set in cross-validation
            classifyResults <- mapply(function(datasetID, classesUse, classifierID, nFeaturesUse) # Could be two or more datasets, classifiers.
            {
              if(datasetID == "dataset")
              {
                measurementsUse <- measurements
              } else { # Use a specific data table.
                if(is(measurements, "DataFrame"))
                  measurementsUse <- measurements[, mcols(measurements[, "dataset"]) == datasetID]
                else # DataFrameList
                  measurementsUse <- measurements[[datasetID]]
              }
              modellingParams <- switch(classifierID, DLDA = ModellingParams(selectParams = SelectParams(tuneParams = list(nFeatures = nFeaturesUse, performanceType = "Balanced Error"))), # Difference in means and DLDA is the default if nothing specified.
                                        FishersLDA = ModellingParams(transformParams = TransformParams(subtractFromLocation, intermediate = "training"), selectParams = SelectParams(leveneRanking, tuneParams = list(nFeatures = nFeaturesUse, performanceType = "Balanced Error")), trainParams = TrainParams(fisherDiscriminant), predictParams = NULL),
                                        NSC = ModellingParams(selectParams = NULL, trainParams = TrainParams(NSCtrainInterface, getFeatures = NSCfeatures), predictParams = PredictParams(NSCpredictInterface)),
                                        naiveBayes = ModellingParams(selectParams = SelectParams(DMDranking, tuneParams = list(nFeatures = nFeaturesUse, performanceType = "Balanced Error")), trainParams = TrainParams(naiveBayesKernel, difference = "weighted"), predictParams = NULL),
                                        elasticNet = ModellingParams(selectParams = NULL, trainParams = TrainParams(elasticNetGLMtrainInterface, getFeatures = elasticNetFeatures), predictParams = PredictParams(elasticNetGLMpredictInterface)),
                                        SVM = ModellingParams(selectParams = NULL, trainParams = TrainParams(SVMtrainInterface, kernel = "polynomial", tuneParams = list(degree = 2:8, cost = 10^(-5:5))), predictParams = PredictParams(SVMpredictInterface)),
                                        randomForest = ModellingParams(selectParams = NULL, trainParams = TrainParams(randomForestTrainInterface, mtry = 0.4 * ncol(measurements)), predictParams = PredictParams(randomForestPredictInterface)),
                                        kTSP = ModellingParams(balancing = "none", selectParams = SelectParams(pairsDifferencesRanking, subsetToSelections = FALSE, featurePairs = interactors, tuneParams = list(nFeatures = nFeaturesUse, performanceType = "Balanced Error")), trainParams = TrainParams(kTSPclassifier, difference = "weighted", intermediate = setNames("selectedFeatures", "featurePairs")), predictParams = NULL) 
              )
              runTests(measurementsUse, classesUse, crossValParams = crossValParams, modellingParams = modellingParams, ...) # Pass through other info such as characteristics data frame.
              
            }, datasetIDs, classes, classifier, nFeatures, SIMPLIFY = FALSE)
            
            if(datasetMode == "allPairs") # Do all cross-classifications.
            { # Test each fitted model on a single data set on all other data sets using the same model.
              extras <- rlang::dots_list(...)
              if("characteristics" %i% names(extras)) characteristicsTable <- extras[["characteristics"]] else characteristicsTable <- DataFrame()
              crossClassifyResults <- mapply(function(trainDatasetID, predictDatasetID, classifier)
              {
                if(predictDatasetID == "dataset")
                {
                  measurementsUse <- measurements
                } else { # Use a specific data table.
                  if(is(measurements, "DataFrame"))
                    measurementsUse <- measurements[, mcols(measurements[, "dataset"]) == predictDatasetID]
                  else # DataFrameList
                    measurementsUse <- measurements[[predictDatasetID]]
                }
                classifyResultUse <- classifyResults[[match(trainDatasetID, datasetIDs)]] # Get existing selections.
                previousSelectParams <- SelectParams(previousSelection, classifyResult = classifyResultUse, intermediate = ".iteration")
                previousSelectNoSubParams <- SelectParams(previousSelection, classifyResult = classifyResultUse, subsetToSelections = FALSE, intermediate = ".iteration")
                modellingParams <- switch(classifierID, DLDA = ModellingParams(selectParams = previousSelectParams), # Difference in means and DLDA is the default if nothing specified.
                                          FishersLDA = ModellingParams(transformParams = TransformParams(subtractFromLocation, intermediate = "training"), selectParams = previousSelectParams, trainParams = TrainParams(fisherDiscriminant), predictParams = NULL),
                                          NSC = ModellingParams(selectParams = previousSelectParams, trainParams = TrainParams(NSCtrainInterface, getFeatures = NSCfeatures), predictParams = PredictParams(NSCpredictInterface)),
                                          naiveBayes = ModellingParams(selectParams = previousSelectParams, trainParams = TrainParams(naiveBayesKernel, difference = "weighted"), predictParams = NULL),
                                          elasticNet = ModellingParams(selectParams = previousSelectParams, trainParams = TrainParams(elasticNetGLMtrainInterface, getFeatures = elasticNetFeatures), predictParams = PredictParams(elasticNetGLMpredictInterface)),
                                          SVM = ModellingParams(selectParams = previousSelectParams, trainParams = TrainParams(SVMtrainInterface, kernel = "polynomial", tuneParams = list(degree = 2:8, cost = 10^(-5:5))), predictParams = PredictParams(SVMpredictInterface)),
                                          randomForest = ModellingParams(selectParams = previousSelectParams, trainParams = TrainParams(randomForestTrainInterface, mtry = 0.4 * ncol(measurements)), predictParams = PredictParams(randomForestPredictInterface)),
                                          kTSP = ModellingParams(balancing = "none", selectParams = previousSelectNoSubParams, trainParams = TrainParams(kTSPclassifier, difference = "weighted", intermediate = setNames("selectedFeatures", "featurePairs")), predictParams = NULL)
              )
              datasetsInfo <- DataFrame(characteristic = c("Training Data Set", "Prediction Data Set"), value = c(trainDatasetID, predictDatasetID))
              characteristicsTable <- rbind(characteristicsTable, datasetsInfo)
              extras[["characteristics"]] <- characteristicsTable
              do.call(runTests, append(list(measurementsUse, classes[predictDatasetID], crossValParams = crossValParams, modellingParams = modellingParams), extras)) # Pass through other info such as characteristics data frame.
              }, diffCombinations[, "trainDatasetID"], diffCombinations[, "predictDatasetID"], diffCombinations[, "classifier"])
              
              # Join individual classifications to cross-classifications.
              classifyResults <- mapply(function(classifyResult, datasetID)
                                 {
                classifyResult@characteristics <- rbind(classifyResult@characteristics, DataFrame(characteristic = c("Training Data Set", "Prediction Data Set"),
                                                                                                  value = rep(datasetIDs, 2)))
                                 })
              classifyResults <- append(classifyResults, crossClassifyResults)
            }
            
            if(length(classifyResults) == 1) # Unpackage it to a ClassifyResult.
              classifyResults <- unlist(classifyResults)
            classifyResults
          })

# One or more omics data sets, possibly with clinical data.
setMethod("quickClassify", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), classes, datasetMode = c("each", "combine", "allPairs"), ...)
          {
            datasetMode <- match.arg(datasetMode)
            tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes, restrict = NULL)
            measurements <- tablesAndClasses[["dataTable"]]
            classes <- tablesAndClasses[["classes"]]
            
            quickClassify(measurements, classes, ...)
          })
