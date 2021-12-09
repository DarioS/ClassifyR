setGeneric("quickClassify", function(measurements, ...)
  standardGeneric("quickClassify"))

setMethod("quickClassify", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
          {
            quickClassify(DataFrame(t(measurements), check.names = FALSE), classes, ...)
          })

setMethod("quickClassify", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes,
crossVal = c("100 permutations, 5 folds", "100 permutations, 25% test", "5-fold", "leave-1-out", "leave-2-out"),
                   seed = 12345, cores = 8, classifier = c("DLDA", "FishersLDA", "NSC", "naiveBayes", "elasticNet", "SVM", "randomForest", "kTSP"),
                  ...)
          {
            if(is.null(rownames(measurements)))
              stop("'measurements' DataFrame must have sample identifiers as its row names.")            
            if(any(is.na(measurements)))
              stop("Some data elements are missing and classifiers don't work with missing data. Consider imputation or filtering.")            
            
            existCores <- parallel::detectCores()
            if(existCores < cores) # Set the number of cores to how many exist.
              cores <- existCores
            if(Sys.info()["sysname"] == "Windows")
            {
              BPparam <- SnowParam(cores, RNGseed = seed)
            } else if (Sys.info()["sysname"] %in% c("MacOS", "Linux")) {
              BPparam <- MulticoreParam(cores, RNGseed = seed) # Multicore is faster than SNOW, but it doesn't work on Windows.
            } else { # Something weird.
              BPparam <- SerialParam()
            }
            
            if(ncol(measurements) < 100)
              nFeatures <- 1:ncol(measurements)
            else if(classifier != "kTSP")
              nFeatures <- seq(10, 100, 10)
            else # k-TSP, so a small number of pairs is desired.
              nFeatures <- 1:10
            crossVal <- match.arg(crossVal)
            
            crossValParams <- switch(crossVal, `100 permutations, 5 folds` = CrossValParams(parallelParams = BPparam, seed = seed),
                             `100 permutations, 25% test` = CrossValParams(samplesSplits = "Permute Percentage Split", parallelParams = BPparam, seed = seed),
                             `5-fold` = CrossValParams(samplesSplits = "k-Fold", parallelParams = BPparam, seed = seed),
                              `leave-2-out` = CrossValParams(samplesSplits = "Leave-k-Out", parallelParams = BPparam, seed = seed),
                              `leave-1-out` = CrossValParams(samplesSplits = "Leave-k-Out", leave = 1, parallelParams = BPparam, seed = seed)
                              )
            classifier <- match.arg(classifier)
            data(HuRI) # Loads a variable called interactors.
            modellingParams <- switch(classifier, DLDA = ModellingParams(selectParams = SelectParams(tuneParams = list(nFeatures = nFeatures, performanceType = "Balanced Error"))), # Difference in means and DLDA is the default if nothing specified.
                                      FishersLDA = ModellingParams(transformParams = TransformParams(subtractFromLocation, intermediate = "training"), selectParams = SelectParams(leveneRanking, tuneParams = list(nFeatures = nFeatures, performanceType = "Balanced Error")), trainParams = TrainParams(fisherDiscriminant), predictParams = NULL),
                                      NSC = ModellingParams(selectParams = NULL, trainParams = TrainParams(NSCtrainInterface, getFeatures = NSCfeatures), predictParams = PredictParams(NSCpredictInterface)),
                                      naiveBayes = ModellingParams(selectParams = SelectParams(DMDrankingtuneParams = list(nFeatures = nFeatures, performanceType = "Balanced Error")), trainParams = TrainParams(naiveBayesKernel, difference = "weighted"), predictParams = NULL),
                                      elasticNet = ModellingParams(selectParams = NULL, trainParams = TrainParams(elasticNetGLMtrainInterface, getFeatures = elasticNetFeatures), predictParams = PredictParams(elasticNetGLMpredictInterface)),
                                      SVM = ModellingParams(selectParams = NULL, trainParams = TrainParams(SVMtrainInterface, kernel = "polynomial", tuneParams = list(degree = 2:8, cost = 10^(-5:5))), predictParams = PredictParams(SVMpredictInterface)),
                                      randomForest = ModellingParams(selectParams = NULL, trainParams = TrainParams(randomForestTrainInterface, mtry = 0.4 * ncol(measurements)), predictParams = PredictParams(randomForestPredictInterface)),
                                      kTSP = ModellingParams(balancing = "none", selectParams = SelectParams(pairsDifferencesRanking, subsetToSelections = FALSE, featurePairs = interactors, tuneParams = list(nFeatures = nFeatures, performanceType = "Balanced Error")), trainParams = TrainParams(kTSPclassifier, difference = "weighted", intermediate = setNames("selectedFeatures", "featurePairs")), predictParams = NULL) 
                                      )
            
            runTests(measurements, classes, crossValParams = crossValParams, modellingParams = modellingParams, ...) # Pass through other info such as characteristics data frame.
          })

# One or more omics data sets, possibly with clinical data.
setMethod("quickClassify", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
          {
            tablesAndClasses <- .MAEtoWideTable(measurements, targets, restrict = NULL)
            measurements <- tablesAndClasses[["dataTable"]]
            classes <- tablesAndClasses[["classes"]]
            
            quickClassify(measurements, classes, ...)
          })