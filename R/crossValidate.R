setGeneric("crossValidate", function(measurements, ...)
  standardGeneric("crossValidate"))

setMethod("crossValidate", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
          {
            crossValidate(DataFrame(t(measurements), check.names = FALSE), classes, ...)
          })

setMethod("crossValidate", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, 
                   classes,
                   nFeatures = NULL,
                   selectionMethod = "t_test",
                   selectionOptimisation = "Resubstitution",
                   classifier = "elasticNet",
                   multiViewMethod = "none",
                   multiViewCombinations = NULL,
                   nFolds = 5,
                   nRepeats = 20, 
                   nCores = 1, 
                   characteristicsLabel = NULL)

          {
            # Check that data is in the right format
            checkData(measurements, 
                      classes)
              
            # Which data-types or data-views are present?
            datasetIDs <- unique(mcols(measurements)[, "dataset"])
           
            # Initiate seed so that comparisons are comparable.
            x <- runif(1)
            seed <- .Random.seed[1]
            
            
            ################################
            #### No multiview
            ################################
            
            if(multiViewMethod == "none"){
            
            # The below loops over dataset and classifier and allows us to answer
            # the following questions:
            #
            # 1) One dataset using one classifier
            # 2) One dataset using multi classifiers
            # 3) Multi datasets individually
            
             resClassifier <- sapply(datasetIDs, function(dataIndex){ # Loop over datasets
               sapply(classifier, function(classifierIndex){ # Loop over classifiers
                 
                 set.seed(seed)
                 CV(measurements = measurements[, mcols(measurements)$dataset == dataIndex], 
                    classes = classes,
                    nFeatures = nFeatures,
                    selectionMethod = selectionMethod,
                    selectionOptimisation = selectionOptimisation,
                    classifier = classifierIndex,
                    multiViewMethod = multiViewMethod,
                    multiViewCombinations = dataIndex,
                    nFolds = nFolds,
                    nRepeats = nRepeats, 
                    nCores = nCores, 
                    characteristicsLabel = characteristicsLabel)
                 },
                 
                 simplify = FALSE)
             },
             
             simplify = FALSE)
             
             result <- unlist(resClassifier)
             
            }
            
            ################################
            #### Yes multiview
            ################################
            
            if(multiViewMethod != "none"){
            if(is.null(multiViewCombinations)) multiViewCombinations <- list(datasetIDs)
            
            dataIndex <- multiViewCombinations[[1]]
            
            result <- CV(measurements = measurements[, mcols(measurements)$dataset %in% dataIndex], 
                         classes = classes,
                         nFeatures = nFeatures,
                         selectionMethod = selectionMethod,
                         selectionOptimisation = selectionOptimisation,
                         classifier = classifier,
                         multiViewMethod = multiViewMethod,
                         multiViewCombinations = dataIndex,
                         nFolds = nFolds,
                         nRepeats = nRepeats, 
                         nCores = nCores, 
                         characteristicsLabel = characteristicsLabel)
            
            }
            
            result
            
          })






# One or more omics data sets, possibly with clinical data.
setMethod("crossValidate", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), classes, datasetMode = c("each", "combine", "allPairs"), ...)
          {
            datasetMode <- match.arg(datasetMode)
            tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes, restrict = NULL)
            measurements <- tablesAndClasses[["dataTable"]]
            classes <- tablesAndClasses[["classes"]]
            
            crossValidate(measurements, classes, ...)
          })






######################################
######################################
generateCrossValParams <- function(nRepeats, nFolds, nCores, selectionOptimisation){
  
  seed <- .Random.seed[1]
  
  if(nCores == 1)
{
    BPparam <- SerialParam(RNGseed = seed)
} else { # Parallel processing is desired.
    # Also set the BPparam RNGseed if the user ran set.seed(someNumber) themselves.
    if(Sys.info()["sysname"] == "Windows") {# Only SnowParam suits Windows.
        BPparam <- BiocParallel::SnowParam(min(nCores, BiocParallel::snowWorkers("SOCK")), RNGseed = seed)
    } else if (Sys.info()["sysname"] %in% c("MacOS", "Linux")) {
        BPparam <- BiocParallel::MulticoreParam(min(nCores, BiocParallel::multicoreWorkers()), RNGseed = seed) # Multicore is faster than SNOW, but it doesn't work on Windows.
    } else { # Something weird.
        BPparam <- BiocParallel::bpparam() # BiocParallel will figure it out.
    }
}
    tuneMode <- selectionOptimisation
    if(tuneMode == "CV") tuneMode <- "Nested CV"
    if(!any(tuneMode %in% c("Resubstitution", "Nested CV"))) stop("selectionOptimisation must be CV or Resubstitution")
    CrossValParams(permutations = nRepeats, folds = nFolds, parallelParams = BPparam, tuneMode = tuneMode)    
}
######################################



######################################
######################################
checkData <- function(measurements, classes){
    if(is.null(rownames(measurements)))
        stop("'measurements' DataFrame must have sample identifiers as its row names.")            
    if(any(is.na(measurements)))
        stop("Some data elements are missing and classifiers don't work with missing data. Consider imputation or filtering.")
    
    # !!!  Need to check mcols has dataset NUm
    
}
######################################



######################################
######################################
generateModellingParams <- function(datasetIDs,
                                    measurements,
                                    nFeatures,
                                    selectionMethod,
                                    selectionOptimisation,
                                    classifier,
                                    multiViewMethod = "none"
                                    ){

    
        
        obsFeatures <- sum(mcols(measurements)[, "dataset"] %in% datasetIDs)
        
        
        if(is.null(nFeatures)) nFeatures <- obsFeatures
        
        if(max(nFeatures) > obsFeatures) {
        
            warning("nFeatures greater than the max number or features in data. 
                                                 Setting to max")
            nFeatures <- pmin(nFeatures, obsFeatures)
        }
        
        classifierName = classifier

        classifier = switch(
            classifier,
            "randomForest" = rfParams(),
            "elasticNet" = elasticParams(),
            "logistic" = logisticParams(),
            "svm" = svmParams(),
            "dlda" = DLDAParams(),
            "naiveBayes" = naiveBayesParams(),
            "elasticNetPreval" = elasticNetPreval()
        )
        
        
        
        selectionMethod = ifelse(is.null(selectionMethod),
                                 "none",
                                 selectionMethod)
        
        selectionMethodParam = switch(
            selectionMethod,
            "none" = ClassifyR::differentMeansRanking,
            "t_test" = ClassifyR::differentMeansRanking,
            "limma" = ClassifyR::limmaSelection,
            "edgeR" = ClassifyR::edgeRselection,
            "NSC" = ClassifyR::NSCselectionInterface,
            "bartlett" = ClassifyR::bartlettSelection,
            "levene" = ClassifyR::leveneSelection,
            "DMD" = ClassifyR::DMDselection,
            "liklihood" = ClassifyR::likelihoodRatioSelection,
            "KS" = ClassifyR::KolmogorovSmirnovSelection,
            "KL" = ClassifyR::KullbackLeiblerSelection
        )
        
        
        selectParams = ClassifyR::SelectParams(
            selectionMethodParam,
            tuneParams = list(nFeatures = nFeatures, 
                              performanceType = "Balanced Error")
            )
        
        params = ModellingParams(
            balancing = "none",
            selectParams = selectParams,
            trainParams = classifier$trainParams,
            predictParams = classifier$predictParams
        )
        
        #if(multiViewMethod != "none") stop("I haven't implemented multiview yet.")
        if(multiViewMethod == "merge"){
            params@trainParams <- TrainParams(mergeTrainInterface, params = params, characteristics = DataFrame(characteristic = "Classifier Name", value = classifierName))
            params@selectParams <- NULL
            }
        # 
        # if(multiViewMethod == "prevalidation"){
        #     params$trainParams <- function(measurements, classes) prevalTrainInterface(measurements, classes, params)
        #     params$trainParams <- function(measurements, classes) prevalTrainInterface(measurements, classes, params)
        # }
        # 
    
            
            params
            
}
######################################






CV <- function(measurements, 
               classes,
               nFeatures = NULL,
               selectionMethod = "t_test",
               selectionOptimisation = "Resubstitution",
               classifier = "elasticNet",
               multiViewMethod = "none",
               multiViewCombinations = NULL,
               nFolds = 5,
               nRepeats = 100, 
               nCores = 1, 
               characteristicsLabel = NULL)

{
  # Check that data is in the right format
  checkData(measurements, 
            classes)
  
  # Which data-types or data-views are present?
  datasetIDs <- unique(mcols(measurements)[, "dataset"])
  if(is.null(multiViewCombinations)) multiViewCombinations <- datasetIDs
  if(is.null(characteristicsLabel)) characteristicsLabel <- "bla"
  
  # Setup cross-validation parameters including 
  crossValParams <- generateCrossValParams(nRepeats = nRepeats, 
                                           nFolds = nFolds, 
                                           nCores = nCores,
                                           selectionOptimisation = selectionOptimisation
  )
  
  # Turn text into TrainParams and TestParams objects
  modellingParams <- generateModellingParams(datasetIDs = datasetIDs,
                                             measurements = measurements,
                                             nFeatures = nFeatures,
                                             selectionMethod = selectionMethod,
                                             selectionOptimisation = selectionOptimisation,
                                             classifier = classifier,
                                             multiViewMethod = multiViewMethod
  )
  
  
  characteristics = S4Vectors::DataFrame(characteristic = c("dataset", "classifier", "selectionMethod", "multiViewMethod", "characteristicsLabel"), value = c(paste(datasetIDs, collapse = ", "), classifier,  selectionMethod, multiViewMethod, characteristicsLabel)) 
  
  classifyResults <- runTests(measurements, classes, crossValParams = crossValParams, modellingParams = modellingParams, characteristics = characteristics)
  classifyResults
  
}











