#'
#' Quick and Easy Cross-validated Classification
#'
#' This convenience function allows specification of cross-validation and
#' classification by a couple of pre-defined keywords and all parameter setting
#' is automatically taken care of. This avoids the need to create any S4
#' parameter objects.
#'
#'
#' @aliases quickClassify quickClassify,matrix-method
#' quickClassify,DataFrameOrDataFrameList-method
#' quickClassify,MultiAssayExperiment-method
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix}, the rows are features, and the columns are samples.
#' @param classes A vector of class labels of class \code{\link{factor}} of the
#' same length as the number of samples in \code{measurements} if it is a
#' \code{\link{matrix}} (i.e. number of columns) or a \code{\link{DataFrame}}
#' (i.e. number of rows) or a character vector of length 1 containing the
#' column name in \code{measurements} if it is a \code{\link{DataFrame}} or the
#' column name in \code{colData(measurements)} if \code{measurements} is a
#' \code{\link{MultiAssayExperiment}}. If a column name, that column will be
#' removed before training.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"clinical"} is also a valid value
#' and specifies that the clinical data table will be used.
#' @param permutations Default: 100. Number of times to permute the samples
#' before folding them. Set to \code{NULL} for a singlee k-fold
#' cross-validation.
#' @param folds Default: 5. The number of disjoint sets to partition the
#' samples into. Set this value to the number of features in
#' \code{measurements} (if it is not a list of data sets of unequal sample
#' numbers) to effectively perform a traditional leave-one-out
#' cross-validation. This variable must be some integer.
#' @param cores Default: 1. Number of CPU cores to use. Each cross-validation
#' iteration is done on one core.
#' @param classifier A short name describing the kind of classification to do.
#' A vector means that each classification will be applied to the data in turn.
#' Specify one name if only one classifier desired. The options represent:
#' \itemize{\item \code{"DLDA"}: Feature selection based on differences in
#' means ranking and diagonal LDA classifier.  \item \code{"FishersLDA"}:
#' Transformation by subtraction of values from the mean value of all samples
#' in the training set followed by selection by Levene's test statistic ranking
#' based on differential variability and classification by Fisher's LDA.  \item
#' \code{"NSC"}: Nearest shrunken centroid classifier. No explicit feature
#' selection before classifier.  \item \code{"naiveBayes"}: Feature selection
#' based on differences in median and Qn followed by naive Bayes classifier
#' using weighted voting based on the difference between class densities. See
#' \code{\link{DMDranking}} for more details. Good for differential
#' dsitribution classification.  \item \code{"elasticNet"}: Elastic Net GLM
#' classifier with \code{"multinomial"} as the family.  \item \code{"SVM"}:
#' Polynomial kernel support vector machine. Does parameter tuning of poly
#' nomial degree between 2 and 8 and the cost penality between 0.00001 to
#' 100000 in power increments of 1.  \item \code{"randomForest"}: Random forest
#' with the number of variables tried at each split eing 40\% of the total
#' number of variables in the data set.  \item \code{"kTSP"}: Feature selection
#' based on pairs differences ranking and classification based on weighted
#' voting of the selected feature pairs.  } \item \code{datasetMode} Default: "each".
#' Whether to classify each input data set individually or concatenate
#' different ones into a single table by specifying \code{"combine"} or do
#' classification on all possible pairs of data sets, using the first data set
#' of the pair for training and the second of the pair for making predictions
#' with by specifying \code{"allPairs"}. Only relevant if \code{measurements}
#' is a \code{DataFrame} created by the \code{MultiAssayExperiment} method,
#' possibly containing variables from multiple assays or a
#' \code{DataFrameList}.  \item...Either variables not used by the
#' \code{matrix} nor the \code{MultiAssayExperiment} method which are passed
#' into and used by the \code{DataFrame} method or variables not used by
#' \code{quickClassify} but used by runTests (i.e. \code{characteristics} and
#' \code{verbose}.
#'
#' @importFrom rlang dots_list
#' @export

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
