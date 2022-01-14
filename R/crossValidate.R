setGeneric("crossValidate", function(measurements, ...)
    standardGeneric("crossValidate"))

setMethod("crossValidate", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
          {
              crossValidate(S4Vectors::DataFrame(t(measurements), check.names = FALSE), classes, ...)
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
              # Check that other variables are in the right format and fix
              nFeatures <- cleanNFeatures(nFeatures = nFeatures,
                                          measurements = measurements)
              selectionMethod <- cleanSelectionMethod(selectionMethod = selectionMethod,
                                                      measurements = measurements)
              classifier <- cleanClassifier(classifier = classifier,
                                            measurements = measurements)
              
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
                  
                  resClassifier <-
                      sapply(datasetIDs, function(dataIndex) {
                          # Loop over datasets
                          sapply(classifier[[dataIndex]], function(classifierIndex) {
                              # Loop over classifiers
                              sapply(selectionMethod[[dataIndex]], function(selectionIndex) {
                                  # Loop over classifiers
                                  
                                  set.seed(seed)
                                  CV(
                                      measurements = measurements[, mcols(measurements)$dataset == dataIndex],
                                      classes = classes,
                                      nFeatures = nFeatures[dataIndex],
                                      selectionMethod = selectionIndex,
                                      selectionOptimisation = selectionOptimisation,
                                      classifier = classifierIndex,
                                      multiViewMethod = multiViewMethod,
                                      multiViewCombinations = dataIndex,
                                      nFolds = nFolds,
                                      nRepeats = nRepeats,
                                      nCores = nCores,
                                      characteristicsLabel = characteristicsLabel
                                  )
                              },
                              
                              simplify = FALSE)
                          },
                          
                          simplify = FALSE)
                      },
                      
                      simplify = FALSE)
                  
                  result <- unlist(unlist(resClassifier))
                  
              }
              
              ################################
              #### Yes multiview
              ################################
              
              ### Merging or binding to combine data
              if(multiViewMethod == "merge"){
                  
                  
                  # The below loops over different combinations of datasets and merges them together.
                  # This allows someone to answer which combinations of the datasets might be most useful.
                  
                  
                  if(is.null(multiViewCombinations)) multiViewCombinations <- do.call("c", sapply(seq_len(length(datasetIDs)),function(n)combn(datasetIDs, n, simplify = FALSE)))
                  
                  
                  result <- sapply(multiViewCombinations, function(dataIndex){
                      CV(measurements = measurements[, mcols(measurements)$dataset %in% dataIndex], 
                         classes = classes,
                         nFeatures = nFeatures[dataIndex],
                         selectionMethod = selectionMethod[dataIndex],
                         selectionOptimisation = selectionOptimisation,
                         classifier = classifier[dataIndex],
                         multiViewMethod = ifelse(length(dataIndex)==1, "none", multiViewMethod),
                         multiViewCombinations = dataIndex,
                         nFolds = nFolds,
                         nRepeats = nRepeats, 
                         nCores = nCores, 
                         characteristicsLabel = characteristicsLabel)
                  }, simplify = FALSE)
                  
              }
              
              
              ### Prevalidation to combine data
              if(multiViewMethod == "prevalidation"){
                  
                  
                  # The below loops over different combinations of datasets and combines them together using prevalidation.
                  # This allows someone to answer which combinations of the datasets might be most useful.
                  
                  
                  if(is.null(multiViewCombinations)){
                      multiViewCombinations <- do.call("c", sapply(seq_len(length(datasetIDs)),function(n)combn(datasetIDs, n, simplify = FALSE)))
                      multiViewCombinations <- multiViewCombinations[sapply(multiViewCombinations, function(x)"clinical"%in%x, simplify = TRUE)]
                      if(length(multiViewCombinations)==0) stop("No multiViewCombinations with `clinical` data")
                  }
                  
                  
                  result <- sapply(multiViewCombinations, function(dataIndex){
                      CV(measurements = measurements[, mcols(measurements)$dataset %in% dataIndex], 
                         classes = classes,
                         nFeatures = nFeatures[dataIndex],
                         selectionMethod = selectionMethod[dataIndex],
                         selectionOptimisation = selectionOptimisation,
                         classifier = classifier[dataIndex],
                         multiViewMethod = ifelse(length(dataIndex)==1, "none", multiViewMethod),
                         multiViewCombinations = dataIndex,
                         nFolds = nFolds,
                         nRepeats = nRepeats, 
                         nCores = nCores, 
                         characteristicsLabel = characteristicsLabel)
                  }, simplify = FALSE)
                  
              }
              
              
              
              ### Prevalidation to combine data
              if(multiViewMethod == "pca"){
                  
                  
                  # The below loops over different combinations of datasets and combines them together using prevalidation.
                  # This allows someone to answer which combinations of the datasets might be most useful.
                  
                  
                  if(is.null(multiViewCombinations)){
                      multiViewCombinations <- do.call("c", sapply(seq_len(length(datasetIDs)),function(n)combn(datasetIDs, n, simplify = FALSE)))
                      multiViewCombinations <- multiViewCombinations[sapply(multiViewCombinations, function(x)"clinical"%in%x, simplify = TRUE)]
                      if(length(multiViewCombinations)==0) stop("No multiViewCombinations with `clinical` data")
                  }
                  
                  
                  result <- sapply(multiViewCombinations, function(dataIndex){
                      CV(measurements = measurements[, mcols(measurements)$dataset %in% dataIndex], 
                         classes = classes,
                         nFeatures = nFeatures[dataIndex],
                         selectionMethod = selectionMethod[dataIndex],
                         selectionOptimisation = selectionOptimisation,
                         classifier = classifier[dataIndex],
                         multiViewMethod = ifelse(length(dataIndex)==1, "none", multiViewMethod),
                         multiViewCombinations = dataIndex,
                         nFolds = nFolds,
                         nRepeats = nRepeats, 
                         nCores = nCores, 
                         characteristicsLabel = characteristicsLabel)
                  }, simplify = FALSE)
                  
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
cleanNFeatures <- function(nFeatures, measurements){
    #### Clean up 
    obsFeatures <- unlist(as.list(table(mcols(measurements)[, "dataset"])))
    if(is.null(names(nFeatures)) & length(nFeatures) == 1) nFeatures <- as.list(pmin(obsFeatures, nFeatures))
    if(is.null(names(nFeatures)) & length(nFeatures) > 1) nFeatures <- sapply(obsFeatures, function(x)pmin(obsFeatures, nFeatures), simplify = FALSE)
    #if(is.null(names(nFeatures)) & length(nFeatures) > 1) stop("nFeatures needs to be a named numeric vector or list with the same names as the datasets.")
    if(!all(names(obsFeatures) %in% names(nFeatures))) stop("nFeatures needs to be a named numeric vector or list with the same names as the datasets.")
    if(all(names(obsFeatures) %in% names(nFeatures)) & is(nFeatures, "numeric")) nFeatures <- as.list(pmin(obsFeatures, nFeatures[names(obsFeatures)]))
    if(all(names(obsFeatures) %in% names(nFeatures)) & is(nFeatures, "list")) nFeatures <- mapply(pmin, nFeatures[names(obsFeatures)], obsFeatures, SIMPLIFY = FALSE)
    nFeatures
}

######################################
######################################
cleanSelectionMethod <- function(selectionMethod, measurements){
    #### Clean up 
    obsFeatures <- unlist(as.list(table(mcols(measurements)[, "dataset"])))

    if(is.null(names(selectionMethod)) & length(selectionMethod) == 1) selectionMethod <- sapply(names(obsFeatures), function(x)selectionMethod, simplify = FALSE)
    if(is.null(names(selectionMethod)) & length(selectionMethod) > 1) selectionMethod <- sapply(names(obsFeatures), function(x)selectionMethod, simplify = FALSE)
    #if(is.null(names(selectionMethod)) & length(selectionMethod) > 1) stop("selectionMethod needs to be a named character vector or list with the same names as the datasets.")
    if(!all(names(obsFeatures) %in% names(selectionMethod))) stop("selectionMethod needs to be a named character vector or list with the same names as the datasets.")
    if(all(names(obsFeatures) %in% names(selectionMethod)) & is(selectionMethod, "character")) selectionMethod <- as.list(selectionMethod[names(obsFeatures)])
    selectionMethod
    
}

######################################
######################################
cleanClassifier <- function(classifier, measurements){
    #### Clean up 
    obsFeatures <- unlist(as.list(table(mcols(measurements)[, "dataset"])))

    if(is.null(names(classifier)) & length(classifier) == 1) classifier <- sapply(names(obsFeatures), function(x)classifier, simplify = FALSE)
    if(is.null(names(classifier)) & length(classifier) > 1) classifier <- sapply(names(obsFeatures), function(x)classifier, simplify = FALSE)
    #if(is.null(names(classifier)) & length(classifier) > 1) stop("classifier needs to be a named character vector or list with the same names as the datasets.")
    if(!all(names(obsFeatures) %in% names(classifier))) stop("classifier needs to be a named character vector or list with the same names as the datasets.")
    if(all(names(obsFeatures) %in% names(classifier)) & is(classifier, "character")) classifier <- as.list(classifier[names(obsFeatures)])
    classifier
    
}


######################################
######################################
#' Title
#'
#' @param nRepeats 
#' @param nFolds 
#' @param nCores 
#' @param selectionOptimisation 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' @import BiocParallel
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
    if(!any(tuneMode %in% c("Resubstitution", "Nested CV", "none"))) stop("selectionOptimisation must be CV or Resubstitution or none")
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
    
    
    if(multiViewMethod != "none") {
        params <- generateMultiviewParams(datasetIDs,
                                          measurements,
                                          nFeatures,
                                          selectionMethod,
                                          selectionOptimisation,
                                          classifier,
                                          multiViewMethod)
        return(params)
    }
    
    
    
    
    obsFeatures <- sum(mcols(measurements)[, "dataset"] %in% datasetIDs)
    
    
    nFeatures <- unlist(nFeatures)
    
    if(max(nFeatures) > obsFeatures) {
        
        warning("nFeatures greater than the max number of features in data. 
                                                 Setting to max")
        nFeatures <- pmin(nFeatures, obsFeatures)
    }
    
    classifier <- unlist(classifier)
    
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
    
    
    selectionMethod <- unlist(selectionMethod)
    
    selectionMethod <- ifelse(is.null(selectionMethod),
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
    
    # 
    # if(multiViewMethod == "prevalidation"){
    #     params$trainParams <- function(measurements, classes) prevalTrainInterface(measurements, classes, params)
    #     params$trainParams <- function(measurements, classes) prevalTrainInterface(measurements, classes, params)
    # }
    # 
    
    
    params
    
}
######################################



#' Title
#'
#' @param datasetIDs 
#' @param measurements 
#' @param nFeatures 
#' @param selectionMethod 
#' @param selectionOptimisation 
#' @param classifier 
#' @param multiViewMethod 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' @importFrom S4Vectors DataFrame 
generateMultiviewParams <- function(datasetIDs,
                                    measurements,
                                    nFeatures,
                                    selectionMethod,
                                    selectionOptimisation,
                                    classifier,
                                    multiViewMethod){
    
    if(multiViewMethod == "merge"){
        
        if(length(classifier) > 1) classifier <- classifier[[1]]
        
        # Split measurements up by dataset.
        assayTrain <- sapply(datasetIDs, function(x) measurements[,mcols(measurements)[["dataset"]]%in%x], simplify = FALSE)
        
        # Generate params for each dataset. This could be extended to have different selectionMethods for each type
        paramsDatasets <- mapply(generateModellingParams,
                                 nFeatures = nFeatures[datasetIDs], 
                                 selectionMethod = selectionMethod[datasetIDs],
                                 datasetIDs = datasetIDs,
                                 measurements = assayTrain[datasetIDs],
                                 MoreArgs = list(
                                     selectionOptimisation = selectionOptimisation,
                                     classifier = classifier,
                                     multiViewMethod = "none"),
                                 SIMPLIFY = FALSE)
        
        # Generate some params for merged model.
        params <- generateModellingParams(datasetIDs = datasetIDs,
                                          measurements = measurements,
                                          nFeatures = nFeatures,
                                          selectionMethod = selectionMethod,
                                          selectionOptimisation = "none",
                                          classifier = classifier,
                                          multiViewMethod = "none")
        
        # Update selectParams to use 
        params@selectParams <- SelectParams(selectMulti, 
                                            params = paramsDatasets, 
                                            characteristics = S4Vectors::DataFrame(characteristic = "Selection Name", value = "merge"),
                                            tuneParams = list(nFeatures = nFeatures[[1]], 
                                                              performanceType = "Balanced Error",
                                                              tuneMode = "none")
        )
        return(params)
    }
    
    if(multiViewMethod == "prevalidation"){
        
        # Split measurements up by dataset.
        assayTrain <- sapply(datasetIDs, function(x) measurements[,mcols(measurements)[["dataset"]]%in%x], simplify = FALSE)
        
        # Generate params for each dataset. This could be extended to have different selectionMethods for each type
        paramsDatasets <- mapply(generateModellingParams,
                                 nFeatures = nFeatures[datasetIDs], 
                                 selectionMethod = selectionMethod[datasetIDs],
                                 datasetIDs = datasetIDs,
                                 measurements = assayTrain[datasetIDs],
                                 classifier = classifier[datasetIDs],
                                 MoreArgs = list(
                                     selectionOptimisation = selectionOptimisation,
                                     multiViewMethod = "none"),
                                 SIMPLIFY = FALSE)
        
        
        params <- ModellingParams(
            balancing = "none",
            selectParams = NULL,
            trainParams = TrainParams(prevalTrainInterface, params = paramsDatasets, characteristics = paramsDatasets$clinical@trainParams@characteristics),
            predictParams = PredictParams(prevalPredictInterface, characteristics = paramsDatasets$clinical@predictParams@characteristics)
        )
        
        return(params)
    }
    
    if(multiViewMethod == "prevalidation"){
        
        # Split measurements up by dataset.
        assayTrain <- sapply(datasetIDs, function(x) measurements[,mcols(measurements)[["dataset"]]%in%x], simplify = FALSE)
        
        # Generate params for each dataset. This could be extended to have different selectionMethods for each type
        paramsDatasets <- mapply(generateModellingParams,
                                 nFeatures = nFeatures[datasetIDs], 
                                 selectionMethod = selectionMethod[datasetIDs],
                                 datasetIDs = datasetIDs,
                                 measurements = assayTrain[datasetIDs],
                                 classifier = classifier[datasetIDs],
                                 MoreArgs = list(
                                     selectionOptimisation = selectionOptimisation,
                                     multiViewMethod = "none"),
                                 SIMPLIFY = FALSE)
        
        
        params <- ModellingParams(
            balancing = "none",
            selectParams = NULL,
            trainParams = TrainParams(prevalTrainInterface, params = paramsDatasets, characteristics = paramsDatasets$clinical@trainParams@characteristics),
            predictParams = PredictParams(prevalPredictInterface, characteristics = paramsDatasets$clinical@predictParams@characteristics)
        )
        
        return(params)
    }
    
    
    if(multiViewMethod == "pca"){
        
        # Split measurements up by dataset.
        assayTrain <- sapply(datasetIDs, function(x) measurements[,mcols(measurements)[["dataset"]]%in%x], simplify = FALSE)
        
        # Generate params for each dataset. This could be extended to have different selectionMethods for each type
        paramsClinical <-  list(clinical = generateModellingParams(
                                 nFeatures = nFeatures["clinical"], 
                                 selectionMethod = selectionMethod["clinical"],
                                 datasetIDs = "clinical",
                                 measurements = assayTrain[["clinical"]],
                                 classifier = classifier["clinical"],
                                 selectionOptimisation = selectionOptimisation,
                                 multiViewMethod = "none"))
        
        
        params <- ModellingParams(
            balancing = "none",
            selectParams = NULL,
            trainParams = TrainParams(pcaTrainInterface, params = paramsClinical, nFeatures = nFeatures, characteristics = paramsClinical$clinical@trainParams@characteristics),
            predictParams = PredictParams(pcaPredictInterface, characteristics = paramsClinical$clinical@predictParams@characteristics)
        )
        
        return(params)
    }
    
}


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
    
    # Check that other variables are in the right format and fix
    nFeatures <- cleanNFeatures(nFeatures = nFeatures,
                                measurements = measurements)
    selectionMethod <- cleanSelectionMethod(selectionMethod = selectionMethod,
                                            measurements = measurements)
    classifier <- cleanClassifier(classifier = classifier,
                                  measurements = measurements)
    
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
    
    
    characteristics = S4Vectors::DataFrame(characteristic = c("dataset", "classifier", "selectionMethod", "multiViewMethod", "characteristicsLabel"), value = c(paste(datasetIDs, collapse = ", "), paste(classifier, collapse = ", "),  paste(selectionMethod, collapse = ", "), multiViewMethod, characteristicsLabel)) 
    
    classifyResults <- runTests(measurements, classes, crossValParams = crossValParams, modellingParams = modellingParams, characteristics = characteristics)
    classifyResults
    
}






simplifyResults <- function(results, values = c("dataset", "classifier", "selectionMethod", "multiViewMethod")){
    ch <- sapply(results, function(x) x@characteristics[x@characteristics$characteristic %in% values, "value"], simplify = TRUE)
    ch <- data.frame(t(ch))
    results[!duplicated(ch)]
}




