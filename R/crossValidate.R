#' Cross-validation to evaluate classification performance.
#' 
#' This function has been designed to faciliate the comparison of classification
#'  methods using cross-validation. A selection of typical 
#' comparisons are implemented.
#'
#' @param measurements Either a \code{\link{DataFrame}}, \code{\link{data.frame}}, \code{\link{matrix}}, \code{\link{MultiAssayExperiment}} 
#' or a list of these objects containing the training data.  For a
#' \code{matrix} and \code{data.frame}, the rows are samples and the columns are features. For a \code{data.frame} or \code{\link{MultiAssayExperiment}}
#' the rows are features and the columns are samples as is typical in Bioconductor.
#' @param classes A vector of class labels of class \code{\link{factor}} of the
#' same length as the number of samples in \code{measurements} or a character vector of length 1 containing the
#' column name in \code{measurements} if it is a \code{\link{DataFrame}} or the
#' column name in \code{colData(measurements)} if \code{measurements} is a \code{\link{MultiAssayExperiment}}. If a column name, that column will be
#' removed before training.
#' @param nFeatures The number of features to be used for classification. If this is a single number, the same number of features will be used for all comparisons
#' or datasets. If a numeric vector these will be optimised over using \code{selectionOptimisation}. If a named vector with the same names of multiple datasets, 
#' a different number of features will be used for each dataset. If a named list of vectors, the respective number of features will be optimised over. 
#' Set to NULL or "all" if all features should be used.
#' @param selectionMethod A character vector of feature selection methods to compare. If a named character vector with names corresponding to different datasets, 
#' and performing multiview classification, the respective classification methods will be used on each dataset.
#' @param selectionOptimisation A character of "Resubstitution", "Nested CV" or "none" specifying the approach used to optimise nFeatures. 
#' @param classifier A character vector of classification methods to compare. If a named character vector with names corresponding to different datasets, 
#' and performing multiview classification, the respective classification methods will be used on each dataset.
#' @param multiViewMethod A character vector specifying the multiview method or data integration approach to use.
#' @param dataCombinations A character vector or list of character vectors proposing the datasets or, in the case of a list, combination of datasets to use
#' with each element being a vector of datasets to combine.
#' @param nFolds A numeric specifying the number of folds to use for cross-validation.
#' @param nRepeats A numeric specifying the the number of repeats or permutations to use for cross-validation.
#' @param nCores A numeric specifying the number of cores used if the user wants to use parallelisation. 
#' @param characteristicsLabel A character specifying an additional label for the cross-validation run.
#'
#' @details
#' \code{selectionMethod} can be any of the following implemented approaches - randomForest, elasticNet, logistic, svm, dlda or naiveBayes. 
#' 
#' \code{classifier} can be any of the following implemented approaches -  none, t_test, limma, edgeR, NSC, bartlette, levene, DMD, likelihood, KS or KL.
#' 
#' \code{multiViewMethod} can take a few different values. Using \code{merge} will merge or bind the datasets after feature selection. 
#'  Using \code{prevlidation} will build prevalidated vectors on all the datasets except the clinical data. There must be a dataset called clinical.
#'  Using \code{pca} will perform pca on each dataset and then merge the top few components with the clinical data. There must be a dataset called clinical.
#'
#' @return An object of class \code{\link{ClassifyResult}}
#' @export
#' @aliases crossValidate crossValidate,matrix-method crossValidate,DataFrame-method
#' crossValidate,MultiAssayExperiment-method, crossValidate,data.frame-method
#' @rdname crossValidate
#' @name crossValidate
#'
#' @examples
#' 
#' data(asthma)
#' 
#' # Compare randomForest and svm classifiers.
#' result <- crossValidate(measurements, classes, classifier = c("randomForest", "svm"))
#' Boxplot(result)
#' 
#' 
#' # Compare performance of different datasets. 
#' # First make a toy example dataset with multiple data types. We'll randomly assign different features to be clinical, gene or protein.
#' set.seed(51773)
#' measurements <- DataFrame(t(measurements))
#' mcols(measurements)$dataset <- c(rep("clinical",20),sample(c("gene", "protein"), ncol(measurements)-20, replace = TRUE))
#' mcols(measurements)$feature <- colnames(measurements)
#' 
#' # We'll use different nFeatures for each dataset. We'll also use repeated cross-validation with 5 repeats for speed in the example.
#' set.seed(51773)
#' result <- crossValidate(measurements, classes, nFeatures = c(clinical = 5, gene = 20, protein = 30), classifier = "randomForest", nRepeats = 5)
#' Boxplot(result)
#' 
#' # Merge different datasets. But we will only do this for two combinations. If dataCombinations is not specified it would attempt all combinations.
#' set.seed(51773)
#' result <- crossValidate(measurements, classes, dataCombinations = list(c("clinical", "protein"), c("clinical", "gene")), multiViewMethod = "merge", nRepeats = 5)
#' Boxplot(resultMerge)
#' 
#' 
#' Boxplot(c(result, resultMerge))
#' 
setGeneric("crossValidate", function(measurements,
                                     classes,
                                     nFeatures = 20,
                                     selectionMethod = "t_test",
                                     selectionOptimisation = "Resubstitution",
                                     classifier = "randomForest",
                                     multiViewMethod = "none",
                                     dataCombinations = NULL,
                                     nFolds = 5,
                                     nRepeats = 20,
                                     nCores = 1,
                                     characteristicsLabel = NULL,
                                     ...)
    standardGeneric("crossValidate"))


setMethod("crossValidate", "DataFrame", 
          function(measurements,
                   classes, 
                   nFeatures = 20,
                   selectionMethod = "t_test",
                   selectionOptimisation = "Resubstitution",
                   classifier = "randomForest",
                   multiViewMethod = "none",
                   dataCombinations = NULL,
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
              
              ##!!!!! Do something with data combinations

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
                  
                  # We should probably transition this to use grid instead.

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
                                      dataCombinations = dataIndex,
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


                  if(is.null(dataCombinations)) dataCombinations <- do.call("c", sapply(seq_len(length(datasetIDs)),function(n)combn(datasetIDs, n, simplify = FALSE)))


                  result <- sapply(dataCombinations, function(dataIndex){
                      CV(measurements = measurements[, mcols(measurements)$dataset %in% dataIndex],
                         classes = classes,
                         nFeatures = nFeatures[dataIndex],
                         selectionMethod = selectionMethod[dataIndex],
                         selectionOptimisation = selectionOptimisation,
                         classifier = classifier[dataIndex],
                         multiViewMethod = ifelse(length(dataIndex)==1, "none", multiViewMethod),
                         dataCombinations = dataIndex,
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


                  if(is.null(dataCombinations)){
                      dataCombinations <- do.call("c", sapply(seq_len(length(datasetIDs)),function(n)combn(datasetIDs, n, simplify = FALSE)))
                      dataCombinations <- dataCombinations[sapply(dataCombinations, function(x)"clinical"%in%x, simplify = TRUE)]
                      if(length(dataCombinations)==0) stop("No dataCombinations with `clinical` data")
                  }


                  result <- sapply(dataCombinations, function(dataIndex){
                      CV(measurements = measurements[, mcols(measurements)$dataset %in% dataIndex],
                         classes = classes,
                         nFeatures = nFeatures[dataIndex],
                         selectionMethod = selectionMethod[dataIndex],
                         selectionOptimisation = selectionOptimisation,
                         classifier = classifier[dataIndex],
                         multiViewMethod = ifelse(length(dataIndex)==1, "none", multiViewMethod),
                         dataCombinations = dataIndex,
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


                  if(is.null(dataCombinations)){
                      dataCombinations <- do.call("c", sapply(seq_len(length(datasetIDs)),function(n)combn(datasetIDs, n, simplify = FALSE)))
                      dataCombinations <- dataCombinations[sapply(dataCombinations, function(x)"clinical"%in%x, simplify = TRUE)]
                      if(length(dataCombinations)==0) stop("No dataCombinations with `clinical` data")
                  }


                  result <- sapply(dataCombinations, function(dataIndex){
                      CV(measurements = measurements[, mcols(measurements)$dataset %in% dataIndex],
                         classes = classes,
                         nFeatures = nFeatures[dataIndex],
                         selectionMethod = selectionMethod[dataIndex],
                         selectionOptimisation = selectionOptimisation,
                         classifier = classifier[dataIndex],
                         multiViewMethod = ifelse(length(dataIndex)==1, "none", multiViewMethod),
                         dataCombinations = dataIndex,
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
          function(measurements,
                   classes, 
                   nFeatures = 20,
                   selectionMethod = "t_test",
                   selectionOptimisation = "Resubstitution",
                   classifier = "randomForest",
                   multiViewMethod = "none",
                   dataCombinations = NULL,
                   nFolds = 5,
                   nRepeats = 20,
                   nCores = 1,
                   characteristicsLabel = NULL)
          {
              targets = names(measurements)
              tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes, restrict = NULL)
              measurements <- tablesAndClasses[["dataTable"]]
              classes <- tablesAndClasses[["classes"]]

              crossValidate(measurements = measurements,
                            classes = classes, 
                            nFeatures = nFeatures,
                            selectionMethod = selectionMethod,
                            selectionOptimisation = selectionOptimisation,
                            classifier = classifier,
                            multiViewMethod = multiViewMethod,
                            dataCombinations = dataCombinations,
                            nFolds = nFolds,
                            nRepeats = nRepeats,
                            nCores = nCores,
                            characteristicsLabel = characteristicsLabel)
          })


setMethod("crossValidate", "data.frame", # data.frame of numeric measurements.
          function(measurements,
                   classes, 
                   nFeatures = 20,
                   selectionMethod = "t_test",
                   selectionOptimisation = "Resubstitution",
                   classifier = "randomForest",
                   multiViewMethod = "none",
                   dataCombinations = NULL,
                   nFolds = 5,
                   nRepeats = 20,
                   nCores = 1,
                   characteristicsLabel = NULL)
          {
              measurements <- S4Vectors::DataFrame(t(measurements), check.names = FALSE)
              mcols(measurements)$dataset <- "dataset"
              mcols(measurements)$feature <- colnames(measurements)
              crossValidate(measurements = measurements,
                            classes = classes, 
                            nFeatures = nFeatures,
                            selectionMethod = selectionMethod,
                            selectionOptimisation = selectionOptimisation,
                            classifier = classifier,
                            multiViewMethod = multiViewMethod,
                            dataCombinations = dataCombinations,
                            nFolds = nFolds,
                            nRepeats = nRepeats,
                            nCores = nCores,
                            characteristicsLabel = characteristicsLabel)
          })

setMethod("crossValidate", "matrix", # Matrix of numeric measurements.
          function(measurements,
                   classes, 
                   nFeatures = 20,
                   selectionMethod = "t_test",
                   selectionOptimisation = "Resubstitution",
                   classifier = "randomForest",
                   multiViewMethod = "none",
                   dataCombinations = NULL,
                   nFolds = 5,
                   nRepeats = 20,
                   nCores = 1,
                   characteristicsLabel = NULL)
          {
              measurements <- S4Vectors::DataFrame(t(measurements), check.names = FALSE)
              mcols(measurements)$dataset <- "dataset"
              mcols(measurements)$feature <- colnames(measurements)
              crossValidate(measurements = measurements,
                            classes = classes, 
                            nFeatures = nFeatures,
                            selectionMethod = selectionMethod,
                            selectionOptimisation = selectionOptimisation,
                            classifier = classifier,
                            multiViewMethod = multiViewMethod,
                            dataCombinations = dataCombinations,
                            nFolds = nFolds,
                            nRepeats = nRepeats,
                            nCores = nCores,
                            characteristicsLabel = characteristicsLabel)

          })


###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
setMethod("crossValidate", "list", # data.frame of numeric measurements.
          function(measurements,
                   classes, 
                   nFeatures = 20,
                   selectionMethod = "t_test",
                   selectionOptimisation = "Resubstitution",
                   classifier = "randomForest",
                   multiViewMethod = "none",
                   dataCombinations = NULL,
                   nFolds = 5,
                   nRepeats = 20,
                   nCores = 1,
                   characteristicsLabel = NULL)
          {
              stop("I still need to implement list of datasets")
          })



######################################
######################################
cleanNFeatures <- function(nFeatures, measurements){
    #### Clean up
    obsFeatures <- unlist(as.list(table(mcols(measurements)[, "dataset"])))
    if(is.null(nFeatures) | nFeatures == "all") nFeatures <- as.list(obsFeatures)
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
#' A function to generate a CrossValParams object
#'
#' @inheritParams crossValidate
#'
#' @return CrossValParams object
#' @export
#'
#' @examples
#' CVparams <- generateCrossValParams(nRepeats = 20, nFolds = 5, nCores = 8, selectionOptimisation = "none")
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
#' A function to generate a ModellingParams object
#'
#' @inheritParams crossValidate
#'
#' @return ModellingParams object
#' @export
#'
#' @examples
#' data(asthma)
#' # First make a toy example dataset with multiple data types. We'll randomly assign different features to be clinical, gene or protein.
#' set.seed(51773)
#' measurements <- DataFrame(t(measurements))
#' mcols(measurements)$dataset <- c(rep("clinical",20),sample(c("gene", "protein"), ncol(measurements)-20, replace = TRUE))
#' mcols(measurements)$feature <- colnames(measurements)
#' modellingParams <- generateModellingParams(datasetIDs = c("clinical", "gene", "protein"),
#'                                           measurements = measurements, 
#'                                           nFeatures = list(clinical = 10, gene = 10, protein = 10),
#'                                           selectionMethod = list(clinical = "t_test", gene = "t_test", protein = "t_test"),
#'                                           selectionOptimisation = "none",
#'                                           classifier = "randomForest",
#'                                           multiViewMethod = "merge")
#' @import BiocParallel
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
        "none" = differentMeansRanking,
        "t_test" = differentMeansRanking,
        "limma" = limmaSelection,
        "edgeR" = edgeRselection,
        "NSC" = NSCselectionInterface,
        "bartlett" = bartlettSelection,
        "levene" = leveneSelection,
        "DMD" = DMDselection,
        "liklihood" = likelihoodRatioSelection,
        "KS" = KolmogorovSmirnovSelection,
        "KL" = KullbackLeiblerSelection
    )


    selectParams = SelectParams(
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
               dataCombinations = NULL,
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
    if(is.null(dataCombinations)) dataCombinations <- datasetIDs
    if(is.null(characteristicsLabel)) characteristicsLabel <- "none"

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




