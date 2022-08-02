#' Cross-validation to evaluate classification performance.
#' 
#' This function has been designed to facilitate the comparison of classification
#' methods using cross-validation. A selection of typical  comparisons are implemented.
#'
#' @param measurements Either a \code{\link{DataFrame}}, \code{\link{data.frame}}, \code{\link{matrix}}, \code{\link{MultiAssayExperiment}} 
#' or a list of these objects containing the training data.  For a
#' \code{matrix} and \code{data.frame}, the rows are samples and the columns are features. For a \code{data.frame} or \code{\link{MultiAssayExperiment}} assay
#' the rows are features and the columns are samples, as is typical in Bioconductor.
#' @param outcome A vector of class labels of class \code{\link{factor}} of the
#' same length as the number of samples in \code{measurements} or a character vector of length 1 containing the
#' column name in \code{measurements} if it is a \code{\link{DataFrame}} or the
#' column name in \code{colData(measurements)} if \code{measurements} is a \code{\link{MultiAssayExperiment}}. If a column name, that column will be
#' removed before training. Or a \code{\link{Surv}} object or a character vector of length 2 or 3 specifying the time and event columns in
#' \code{measurements} for survival outcome.
#' @param ... Arguments other than measurements and outcome in the generic.
#' @param nFeatures The number of features to be used for classification. If this is a single number, the same number of features will be used for all comparisons
#' or assays. If a numeric vector these will be optimised over using \code{selectionOptimisation}. If a named vector with the same names of multiple assays, 
#' a different number of features will be used for each assay. If a named list of vectors, the respective number of features will be optimised over. 
#' Set to NULL or "all" if all features should be used.
#' @param selectionMethod A character vector of feature selection methods to compare. If a named character vector with names corresponding to different assays, 
#' and performing multiview classification, the respective classification methods will be used on each assay.
#' @param selectionOptimisation A character of "Resubstitution", "Nested CV" or "none" specifying the approach used to optimise nFeatures. 
#' @param classifier A character vector of classification methods to compare. If a named character vector with names corresponding to different assays, 
#' and performing multiview classification, the respective classification methods will be used on each assay.
#' @param multiViewMethod A character vector specifying the multiview method or data integration approach to use.
#' @param assayCombinations A character vector or list of character vectors proposing the assays or, in the case of a list, combination of assays to use
#' with each element being a vector of assays to combine.
#' @param nFolds A numeric specifying the number of folds to use for cross-validation.
#' @param nRepeats A numeric specifying the the number of repeats or permutations to use for cross-validation.
#' @param nCores A numeric specifying the number of cores used if the user wants to use parallelisation. 
#' @param characteristicsLabel A character specifying an additional label for the cross-validation run.
#' @param object A trained model to predict with.
#' @param newData The data to use to make predictions with.
#'
#' @details
#' \code{classifier} can be any of the following implemented approaches - randomForest, GLM, elasticNetGLM, logistic, SVM, DLDA, kNN, naiveBayes, mixturesNormals.
#' 
#' \code{selectionMethod} can be any of the following implemented approaches -  none, t-test, limma, edgeR, NSC, Bartlett, Levene, DMD, likelihoodRatio, KS or KL.
#' 
#' \code{multiViewMethod} can take a few different values. Using \code{merge} will merge or bind the assays after feature selection. 
#'  Using \code{prevalidation} will build prevalidated vectors on all the assays except the clinical data. There must be a assay called clinical.
#'  Using \code{PCA} will perform Principal Components Analysis on each assay and then merge the top few components with the clinical data. There must be a assay called clinical.
#'
#' @return An object of class \code{\link{ClassifyResult}}
#' @export
#' @aliases crossValidate crossValidate,matrix-method crossValidate,DataFrame-method
#' crossValidate,MultiAssayExperiment-method, crossValidate,data.frame-method
#' @rdname crossValidate
#'
#' @examples
#' 
#' data(asthma)
#' 
#' # Compare randomForest and SVM classifiers.
#' result <- crossValidate(measurements, classes, classifier = c("randomForest", "SVM"))
#' performancePlot(result)
#' 
#' 
#' # Compare performance of different assays. 
#' # First make a toy example assay with multiple data types. We'll randomly assign different features to be clinical, gene or protein.
#' # set.seed(51773)
#' # measurements <- DataFrame(measurements, check.names = FALSE)
#' # mcols(measurements)$assay <- c(rep("clinical",20),sample(c("gene", "protein"), ncol(measurements)-20, replace = TRUE))
#' # mcols(measurements)$feature <- colnames(measurements)
#' 
#' # We'll use different nFeatures for each assay. We'll also use repeated cross-validation with 5 repeats for speed in the example.
#' # set.seed(51773)
#' #result <- crossValidate(measurements, classes, nFeatures = c(clinical = 5, gene = 20, protein = 30), classifier = "randomForest", nRepeats = 5)
#' # performancePlot(result)
#' 
#' # Merge different assays. But we will only do this for two combinations. If assayCombinations is not specified it would attempt all combinations.
#' # set.seed(51773)
#' # resultMerge <- crossValidate(measurements, classes, assayCombinations = list(c("clinical", "protein"), c("clinical", "gene")), multiViewMethod = "merge", nRepeats = 5)
#' # performancePlot(resultMerge)
#' 
#' 
#' # performancePlot(c(result, resultMerge))
#' 
#' @importFrom survival Surv
setGeneric("crossValidate", function(measurements, outcome, ...)
    standardGeneric("crossValidate"))

#' @rdname crossValidate
#' @export
setMethod("crossValidate", "DataFrame", 
          function(measurements,
                   outcome,
                   nFeatures = 20,
                   selectionMethod = "t-test",
                   selectionOptimisation = "Resubstitution",
                   classifier = "randomForest",
                   multiViewMethod = "none",
                   assayCombinations = NULL,
                   nFolds = 5,
                   nRepeats = 20,
                   nCores = 1,
                   characteristicsLabel = NULL)

          {
              # Check that data is in the right format
              splitAssay <- .splitDataAndOutcome(measurements, outcome)
              measurements <- splitAssay[["measurements"]]
              outcome <- splitAssay[["outcome"]]
              
              # Which data-types or data-views are present?
              assayIDs <- unique(mcols(measurements)[, "assay"])
              if(is.null(assayIDs))
                assayIDs <- 1
              
              checkData(measurements, outcome)

              # Check that other variables are in the right format and fix
              nFeatures <- cleanNFeatures(nFeatures = nFeatures,
                                          measurements = measurements)
              selectionMethod <- cleanSelectionMethod(selectionMethod = selectionMethod,
                                                      measurements = measurements)
              if(nFeatures == 1 && classifier == "elasticNetGLM")
              {
                  options(warn = 1)
                  warning("Elastic Net GLM requires two or more features as input but there is only one.
  Using an ordinary GLM instead.")
                  classifier <- "GLM"
              }
              
              classifier <- cleanClassifier(classifier = classifier,
                                            measurements = measurements)
              
              ##!!!!! Do something with data combinations

              # Initiate seed so that comparisons are comparable.
              x <- runif(1)
              seed <- .Random.seed[1]


              ################################
              #### No multiview
              ################################

              if(multiViewMethod == "none"){

                  # The below loops over assay and classifier and allows us to answer
                  # the following questions:
                  #
                  # 1) One assay using one classifier
                  # 2) One assay using multi classifiers
                  # 3) Multi assays individually
                  
                  # We should probably transition this to use grid instead.
                  
                  resClassifier <-
                      sapply(assayIDs, function(assayIndex) {
                          # Loop over assays
                          sapply(classifier[[assayIndex]], function(classifierIndex) {
                              # Loop over classifiers
                              sapply(selectionMethod[[assayIndex]], function(selectionIndex) {
                                  # Loop over classifiers
                                  set.seed(seed)
                                  measurementsUse <- measurements
                                  if(assayIndex != 1) measurementsUse <- measurements[, mcols(measurements)[, "assay"] == assayIndex, drop = FALSE]
                                  CV(
                                      measurements = measurementsUse, outcome = outcome,
                                      assayIDs = assayIndex,
                                      nFeatures = nFeatures[assayIndex],
                                      selectionMethod = selectionIndex,
                                      selectionOptimisation = selectionOptimisation,
                                      classifier = classifierIndex,
                                      multiViewMethod = multiViewMethod,
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


                  # The below loops over different combinations of assays and merges them together.
                  # This allows someone to answer which combinations of the assays might be most useful.


                  if(is.null(assayCombinations)) assayCombinations <- do.call("c", sapply(seq_along(assayIDs), function(nChoose) combn(assayIDs, nChoose, simplify = FALSE)))

                  result <- sapply(assayCombinations, function(assayIndex){
                      CV(measurements = measurements[, mcols(measurements)[["assay"]] %in% assayIndex],
                         outcome = outcome, assayIDs = assayIndex,
                         nFeatures = nFeatures[assayIndex],
                         selectionMethod = selectionMethod[assayIndex],
                         selectionOptimisation = selectionOptimisation,
                         classifier = classifier[assayIndex],
                         multiViewMethod = ifelse(length(assayIndex) == 1, "none", multiViewMethod),
                         nFolds = nFolds,
                         nRepeats = nRepeats,
                         nCores = nCores,
                         characteristicsLabel = characteristicsLabel)
                  }, simplify = FALSE)

              }


              ### Prevalidation to combine data
              if(multiViewMethod == "prevalidation"){


                  # The below loops over different combinations of assays and combines them together using prevalidation.
                  # This allows someone to answer which combinations of the assays might be most useful.


                  if(is.null(assayCombinations))
                  {
                      assayCombinations <- do.call("c", sapply(seq_along(assayIDs), function(nChoose) combn(assayIDs, nChoose, simplify = FALSE)))
                      assayCombinations <- assayCombinations[sapply(assayCombinations, function(combination) "clinical" %in% combination, simplify = TRUE)]
                      if(length(assayCombinations) == 0) stop("No assayCombinations with \"clinical\" data")
                  }


                  result <- sapply(assayCombinations, function(assayIndex){
                      CV(measurements = measurements[, mcols(measurements)[["assay"]] %in% assayIndex],
                         outcome = outcome, assayIDs = assayIndex,
                         nFeatures = nFeatures[assayIndex],
                         selectionMethod = selectionMethod[assayIndex],
                         selectionOptimisation = selectionOptimisation,
                         classifier = classifier[assayIndex],
                         multiViewMethod = ifelse(length(assayIndex) == 1, "none", multiViewMethod),
                         nFolds = nFolds,
                         nRepeats = nRepeats,
                         nCores = nCores,
                         characteristicsLabel = characteristicsLabel)
                  }, simplify = FALSE)

              }



              ### Principal Components Analysis to combine data
              if(multiViewMethod == "PCA"){


                  # The below loops over different combinations of assays and combines them together using prevalidation.
                  # This allows someone to answer which combinations of the assays might be most useful.


                  if(is.null(assayCombinations)){
                      assayCombinations <- do.call("c", sapply(seq_along(assayIDs),function(nChoose) combn(assayIDs, nChoose, simplify = FALSE)))
                      assayCombinations <- assayCombinations[sapply(assayCombinations, function(combination) "clinical" %in% combination, simplify = TRUE)]
                      if(length(assayCombinations) == 0) stop("No assayCombinations with \"clinical\" data")
                  }


                  result <- sapply(assayCombinations, function(assayIndex){
                      CV(measurements = measurements[, mcols(measurements)$assay %in% assayIndex],
                         outcome = outcome, assayIDs = assayIndex,
                         nFeatures = nFeatures[assayIndex],
                         selectionMethod = selectionMethod[assayIndex],
                         selectionOptimisation = selectionOptimisation,
                         classifier = classifier[assayIndex],
                         multiViewMethod = ifelse(length(assayIndex) == 1, "none", multiViewMethod),
                         nFolds = nFolds,
                         nRepeats = nRepeats,
                         nCores = nCores,
                         characteristicsLabel = characteristicsLabel)
                  }, simplify = FALSE)

              }

              result

          })


#' @rdname crossValidate
#' @export
# One or more omics data sets, possibly with clinical data.
setMethod("crossValidate", "MultiAssayExperiment",
          function(measurements,
                   outcome, 
                   nFeatures = 20,
                   selectionMethod = "t-test",
                   selectionOptimisation = "Resubstitution",
                   classifier = "randomForest",
                   multiViewMethod = "none",
                   assayCombinations = NULL,
                   nFolds = 5,
                   nRepeats = 20,
                   nCores = 1,
                   characteristicsLabel = NULL)
          {
              targets <- c(names(measurements), "sampleInfo")
              omicsTargets <- setdiff(targets, "sampleInfo")  
              if(length(omicsTargets) > 0)
              {
                  if(any(anyReplicated(measurements[, , omicsTargets])))
                      stop("Data set contains replicates. Please provide remove or average replicate observations and try again.")
              }
              
              tablesAndoutcome <- .MAEtoWideTable(measurements, targets, outcome, restrict = NULL)
              measurements <- tablesAndoutcome[["dataTable"]]
              outcome <- tablesAndoutcome[["outcome"]]

              crossValidate(measurements = measurements,
                            outcome = outcome, 
                            nFeatures = nFeatures,
                            selectionMethod = selectionMethod,
                            selectionOptimisation = selectionOptimisation,
                            classifier = classifier,
                            multiViewMethod = multiViewMethod,
                            assayCombinations = assayCombinations,
                            nFolds = nFolds,
                            nRepeats = nRepeats,
                            nCores = nCores,
                            characteristicsLabel = characteristicsLabel)
          })

#' @rdname crossValidate
#' @export
setMethod("crossValidate", "data.frame", # data.frame of numeric measurements.
          function(measurements,
                   outcome, 
                   nFeatures = 20,
                   selectionMethod = "t-test",
                   selectionOptimisation = "Resubstitution",
                   classifier = "randomForest",
                   multiViewMethod = "none",
                   assayCombinations = NULL,
                   nFolds = 5,
                   nRepeats = 20,
                   nCores = 1,
                   characteristicsLabel = NULL)
          {
              measurements <- DataFrame(measurements)
              crossValidate(measurements = measurements,
                            outcome = outcome,
                            nFeatures = nFeatures,
                            selectionMethod = selectionMethod,
                            selectionOptimisation = selectionOptimisation,
                            classifier = classifier,
                            multiViewMethod = multiViewMethod,
                            assayCombinations = assayCombinations,
                            nFolds = nFolds,
                            nRepeats = nRepeats,
                            nCores = nCores,
                            characteristicsLabel = characteristicsLabel)
          })

#' @rdname crossValidate
#' @export
setMethod("crossValidate", "matrix", # Matrix of numeric measurements.
          function(measurements,
                   outcome,
                   nFeatures = 20,
                   selectionMethod = "t-test",
                   selectionOptimisation = "Resubstitution",
                   classifier = "randomForest",
                   multiViewMethod = "none",
                   assayCombinations = NULL,
                   nFolds = 5,
                   nRepeats = 20,
                   nCores = 1,
                   characteristicsLabel = NULL)
          {
              measurements <- S4Vectors::DataFrame(measurements, check.names = FALSE)
              crossValidate(measurements = measurements,
                            outcome = outcome,
                            nFeatures = nFeatures,
                            selectionMethod = selectionMethod,
                            selectionOptimisation = selectionOptimisation,
                            classifier = classifier,
                            multiViewMethod = multiViewMethod,
                            assayCombinations = assayCombinations,
                            nFolds = nFolds,
                            nRepeats = nRepeats,
                            nCores = nCores,
                            characteristicsLabel = characteristicsLabel)

          })


###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#' @rdname crossValidate
#' @export
setMethod("crossValidate", "list",
          function(measurements,
                   outcome, 
                   nFeatures = 20,
                   selectionMethod = "t-test",
                   selectionOptimisation = "Resubstitution",
                   classifier = "randomForest",
                   multiViewMethod = "none",
                   assayCombinations = NULL,
                   nFolds = 5,
                   nRepeats = 20,
                   nCores = 1,
                   characteristicsLabel = NULL)
          {
              # Check if the list only contains one data type
              if (measurements |> sapply(class) |> unique() |> length() != 1) {
                  stop("All assays must be of the same type (e.g. data.frame, matrix)")
              }
              
              # Check data type is valid
              if (!(measurements[[1]] |> class() %in% c("data.frame", "DataFrame", "matrix"))) {
                  stop("assays must be of type data.frame, DataFrame or matrix")
              }
              
              # Check the list is named
              if (names(measurements) |> is.null()) {
                  stop("Measurements must be a named list")
              }
              
              # Check same number of samples for all datasets
              if ((measurements |> sapply(dim))[1,] |> unique() |> length() != 1) {
                  stop("All datasets must have the same number of samples")
              }
              
              # Check the number of outcome is the same
              if (((measurements[[1]] |> dim())[1] != length(outcome)) & length(outcome)>2) {
                  stop("outcome must have same number of samples as measurements")
              }
              
              df_list <- sapply(measurements, t, simplify = FALSE)
              df_list <- sapply(measurements , S4Vectors::DataFrame)
              
              df_list <- mapply(function(meas, nam){
                  mcols(meas)$assay <- nam
                  mcols(meas)$feature <- colnames(meas)
                  meas
              }, df_list, names(df_list))
              
              
              combined_df <- do.call(cbind, df_list)
              colnames(combined_df) <- mcols(combined_df)$feature
              
              crossValidate(measurements = combined_df,
                            outcome = outcome, 
                            nFeatures = nFeatures,
                            selectionMethod = selectionMethod,
                            selectionOptimisation = selectionOptimisation,
                            classifier = classifier,
                            multiViewMethod = multiViewMethod,
                            assayCombinations = assayCombinations,
                            nFolds = nFolds,
                            nRepeats = nRepeats,
                            nCores = nCores,
                            characteristicsLabel = characteristicsLabel)
          })



######################################
######################################
cleanNFeatures <- function(nFeatures, measurements){
    #### Clean up
    if(!is.null(mcols(measurements)))
      obsFeatures <- unlist(as.list(table(mcols(measurements)[, "assay"])))
    else obsFeatures <- ncol(measurements)
    if(is.null(nFeatures) || length(nFeatures) == 1 && nFeatures == "all") nFeatures <- as.list(obsFeatures)
    if(is.null(names(nFeatures)) & length(nFeatures) == 1) nFeatures <- as.list(pmin(obsFeatures, nFeatures))
    if(is.null(names(nFeatures)) & length(nFeatures) > 1) nFeatures <- sapply(obsFeatures, function(x)pmin(x, nFeatures), simplify = FALSE)
    #if(is.null(names(nFeatures)) & length(nFeatures) > 1) stop("nFeatures needs to be a named numeric vector or list with the same names as the assays.")
    if(!is.null(names(obsFeatures)) && !all(names(obsFeatures) %in% names(nFeatures))) stop("nFeatures needs to be a named numeric vector or list with the same names as the assays.")
    if(!is.null(names(obsFeatures)) && all(names(obsFeatures) %in% names(nFeatures)) & is(nFeatures, "numeric")) nFeatures <- as.list(pmin(obsFeatures, nFeatures[names(obsFeatures)]))
    if(!is.null(names(obsFeatures)) && all(names(obsFeatures) %in% names(nFeatures)) & is(nFeatures, "list")) nFeatures <- mapply(pmin, nFeatures[names(obsFeatures)], obsFeatures, SIMPLIFY = FALSE)
    nFeatures
}

######################################
######################################
cleanSelectionMethod <- function(selectionMethod, measurements){
    #### Clean up
    if(!is.null(mcols(measurements)))
      obsFeatures <- unlist(as.list(table(mcols(measurements)[, "assay"])))
    else return(list(selectionMethod))

    if(is.null(names(selectionMethod)) & length(selectionMethod) == 1 & !is.null(names(obsFeatures))) selectionMethod <- sapply(names(obsFeatures), function(x) selectionMethod, simplify = FALSE)
    if(is.null(names(selectionMethod)) & length(selectionMethod) > 1 & !is.null(names(obsFeatures))) selectionMethod <- sapply(names(obsFeatures), function(x) selectionMethod, simplify = FALSE)
    #if(is.null(names(selectionMethod)) & length(selectionMethod) > 1) stop("selectionMethod needs to be a named character vector or list with the same names as the assays.")
    if(!is.null(names(obsFeatures)) && !all(names(obsFeatures) %in% names(selectionMethod))) stop("selectionMethod needs to be a named character vector or list with the same names as the assays.")
    if(!is.null(names(obsFeatures)) && all(names(obsFeatures) %in% names(selectionMethod)) & is(selectionMethod, "character")) selectionMethod <- as.list(selectionMethod[names(obsFeatures)])
    selectionMethod
}

######################################
######################################
cleanClassifier <- function(classifier, measurements){
    #### Clean up
    if(!is.null(mcols(measurements)))
      obsFeatures <- unlist(as.list(table(mcols(measurements)[, "assay"])))
    else return(list(classifier))

    if(is.null(names(classifier)) & length(classifier) == 1 & !is.null(names(obsFeatures))) classifier <- sapply(names(obsFeatures), function(x)classifier, simplify = FALSE)
    if(is.null(names(classifier)) & length(classifier) > 1 & !is.null(names(obsFeatures))) classifier <- sapply(names(obsFeatures), function(x)classifier, simplify = FALSE)
    #if(is.null(names(classifier)) & length(classifier) > 1) stop("classifier needs to be a named character vector or list with the same names as the assays.")
    if(!is.null(names(obsFeatures)) && !all(names(obsFeatures) %in% names(classifier))) stop("classifier needs to be a named character vector or list with the same names as the assays.")
    if(!is.null(names(obsFeatures)) && all(names(obsFeatures) %in% names(classifier)) & is(classifier, "character")) classifier <- as.list(classifier[names(obsFeatures)])
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
checkData <- function(measurements, outcome){
    if(is.null(rownames(measurements)))
        stop("'measurements' DataFrame must have sample identifiers as its row names.")
    if(any(is.na(measurements)))
        stop("Some data elements are missing and classifiers don't work with missing data. Consider imputation or filtering.")

    # !!!  Need to check mcols has assay NUm

}
######################################



######################################
######################################
#' A function to generate a ModellingParams object
#'
#' @inheritParams crossValidate
#' @param assayIDs A vector of data set identifiers as long at the number of data sets.
#'
#' @return ModellingParams object
#' @export
#'
#' @examples
#' data(asthma)
#' # First make a toy example assay with multiple data types. We'll randomly assign different features to be clinical, gene or protein.
#' set.seed(51773)
#' measurements <- DataFrame(measurements, check.names = FALSE) 
#' mcols(measurements)$assay <- c(rep("clinical",20),sample(c("gene", "protein"), ncol(measurements)-20, replace = TRUE))
#' mcols(measurements)$feature <- colnames(measurements)
#' modellingParams <- generateModellingParams(assayIDs = c("clinical", "gene", "protein"),
#'                                           measurements = measurements, 
#'                                           nFeatures = list(clinical = 10, gene = 10, protein = 10),
#'                                           selectionMethod = list(clinical = "t-test", gene = "t-test", protein = "t-test"),
#'                                           selectionOptimisation = "none",
#'                                           classifier = "randomForest",
#'                                           multiViewMethod = "merge")
#' @import BiocParallel
generateModellingParams <- function(assayIDs,
                                    measurements,
                                    nFeatures,
                                    selectionMethod,
                                    selectionOptimisation,
                                    classifier,
                                    multiViewMethod = "none"
){
    if(multiViewMethod != "none") {
        params <- generateMultiviewParams(assayIDs,
                                          measurements,
                                          nFeatures,
                                          selectionMethod,
                                          selectionOptimisation,
                                          classifier,
                                          multiViewMethod)
        return(params)
    }




    if(length(assayIDs) > 1) obsFeatures <- sum(mcols(measurements)[, "assay"] %in% assayIDs)
    else obsFeatures <- ncol(measurements)


    nFeatures <- unlist(nFeatures)

    if(max(nFeatures) > obsFeatures) {

        warning("nFeatures greater than the max number of features in data.
                                                 Setting to max")
        nFeatures <- pmin(nFeatures, obsFeatures)
    }

    classifier <- unlist(classifier)

    performanceType <- ifelse(classifier %in% c("CoxPH", "CoxNet", "randomSurvivalForest"), "C-index", "Balanced Accuracy")
    
    
    classifiers <- c("randomForest", "randomSurvivalForest", "GLM", "elasticNetGLM", "SVM", "DLDA",
                     "naiveBayes", "mixturesNormals", "kNN",
                     "CoxPH", "CoxNet")
    # Check classifier
    if(!classifier %in% classifiers)
        stop(paste("Classifier must exactly match of these (be careful of case):", paste(classifiers, collapse = ", ")))
    
    classifier <- switch(
        classifier,
        "randomForest" = RFparams(),
        "randomSurvivalForest" = RSFparams(),
        "GLM" = GLMparams(),
        "elasticNetGLM" = elasticNetGLMparams(),
        "SVM" = SVMparams(),
        "DLDA" = DLDAparams(),
        "naiveBayes" = naiveBayesParams(),
        "mixturesNormals" = mixModelsParams(),
        "kNN" = kNNparams(),
        "CoxPH" = coxphParams(),
        "CoxNet" = coxnetParams()
    )

    selectionMethod <- unlist(selectionMethod)

    selectionMethod <- ifelse(is.null(selectionMethod),
                              "none",
                              selectionMethod)

    selectionMethodParam <- switch(
        selectionMethod,
        "none" = differentMeansRanking,
        "t-test" = differentMeansRanking,
        "limma" = limmaRanking,
        "edgeR" = edgeRranking,
        "Bartlett" = bartlettRanking,
        "Levene" = leveneRanking,
        "DMD" = DMDranking,
        "likelihoodRatio" = likelihoodRatioRanking,
        "KS" = KolmogorovSmirnovRanking,
        "KL" = KullbackLeiblerRanking,
        "CoxPH" = coxphRanking
    )

    selectParams = SelectParams(
        selectionMethodParam,
        tuneParams = list(nFeatures = nFeatures, performanceType = performanceType)
        )
    
    if(selectionMethod == "none" || is.null(selectionMethod)) selectParams <- NULL

    params <- ModellingParams(
        balancing = "none",
        selectParams = selectParams,
        trainParams = classifier$trainParams,
        predictParams = classifier$predictParams
    )

    #if(multiViewMethod != "none") stop("I haven't implemented multiview yet.")

    #
    # if(multiViewMethod == "prevalidation"){
    #     params$trainParams <- function(measurements, outcome) prevalTrainInterface(measurements, outcome, params)
    #     params$trainParams <- function(measurements, outcome) prevalTrainInterface(measurements, outcome, params)
    # }
    #


    params

}
######################################



generateMultiviewParams <- function(assayIDs,
                                    measurements,
                                    nFeatures,
                                    selectionMethod,
                                    selectionOptimisation,
                                    classifier,
                                    multiViewMethod){

    if(multiViewMethod == "merge"){

        if(length(classifier) > 1) classifier <- classifier[[1]]

        # Split measurements up by assay.
        assayTrain <- sapply(assayIDs, function(assayID) if(assayID == 1) measurements else measurements[, mcols(measurements)[["assay"]] %in% assayID], simplify = FALSE)

        # Generate params for each assay. This could be extended to have different selectionMethods for each type
        paramsassays <- mapply(generateModellingParams,
                                 nFeatures = nFeatures[assayIDs],
                                 selectionMethod = selectionMethod[assayIDs],
                                 assayIDs = assayIDs,
                                 measurements = assayTrain[assayIDs],
                                 MoreArgs = list(
                                     selectionOptimisation = selectionOptimisation,
                                     classifier = classifier,
                                     multiViewMethod = "none"),
                                 SIMPLIFY = FALSE)

        # Generate some params for merged model.
        params <- generateModellingParams(assayIDs = assayIDs,
                                          measurements = measurements,
                                          nFeatures = nFeatures,
                                          selectionMethod = selectionMethod,
                                          selectionOptimisation = "none",
                                          classifier = classifier,
                                          multiViewMethod = "none")

        # Update selectParams to use
        params@selectParams <- SelectParams(selectMulti,
                                            params = paramsassays,
                                            characteristics = S4Vectors::DataFrame(characteristic = "Selection Name", value = "merge"),
                                            tuneParams = list(nFeatures = nFeatures[[1]],
                                                              performanceType = "Balanced Error",
                                                              tuneMode = "none")
        )
        return(params)
    }

    if(multiViewMethod == "prevalidation"){

        # Split measurements up by assay.
        assayTrain <- sapply(assayIDs, function(assayID) measurements[, mcols(measurements)[["assay"]] %in% assayID], simplify = FALSE)

        # Generate params for each assay. This could be extended to have different selectionMethods for each type
        paramsassays <- mapply(generateModellingParams,
                                 nFeatures = nFeatures[assayIDs],
                                 selectionMethod = selectionMethod[assayIDs],
                                 assayIDs = assayIDs,
                                 measurements = assayTrain[assayIDs],
                                 classifier = classifier[assayIDs],
                                 MoreArgs = list(
                                     selectionOptimisation = selectionOptimisation,
                                     multiViewMethod = "none"),
                                 SIMPLIFY = FALSE)


        params <- ModellingParams(
            balancing = "none",
            selectParams = NULL,
            trainParams = TrainParams(prevalTrainInterface, params = paramsassays, characteristics = paramsassays$clinical@trainParams@characteristics),
            predictParams = PredictParams(prevalPredictInterface, characteristics = paramsassays$clinical@predictParams@characteristics)
        )

        return(params)
    }

    if(multiViewMethod == "prevalidation"){

        # Split measurements up by assay.
        assayTrain <- sapply(assayIDs, function(assayID) measurements[, mcols(measurements)[["assay"]] %in% assayID], simplify = FALSE)

        # Generate params for each assay. This could be extended to have different selectionMethods for each type
        paramsassays <- mapply(generateModellingParams,
                                 nFeatures = nFeatures[assayIDs],
                                 selectionMethod = selectionMethod[assayIDs],
                                 assayIDs = assayIDs,
                                 measurements = assayTrain[assayIDs],
                                 classifier = classifier[assayIDs],
                                 MoreArgs = list(
                                     selectionOptimisation = selectionOptimisation,
                                     multiViewMethod = "none"),
                                 SIMPLIFY = FALSE)


        params <- ModellingParams(
            balancing = "none",
            selectParams = NULL,
            trainParams = TrainParams(prevalTrainInterface, params = paramsassays, characteristics = paramsassays$clinical@trainParams@characteristics),
            predictParams = PredictParams(prevalPredictInterface, characteristics = paramsassays$clinical@predictParams@characteristics)
        )

        return(params)
    }


    if(multiViewMethod == "PCA"){

        # Split measurements up by assay.
        assayTrain <- sapply(assayIDs, function(assayID) measurements[,mcols(measurements)[["assay"]] %in% assayID], simplify = FALSE)

        # Generate params for each assay. This could be extended to have different selectionMethods for each type
        paramsClinical <-  list(clinical = generateModellingParams(
                                 nFeatures = nFeatures["clinical"],
                                 selectionMethod = selectionMethod["clinical"],
                                 assayIDs = "clinical",
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
               outcome,
               assayIDs,
               nFeatures = NULL,
               selectionMethod = "t-test",
               selectionOptimisation = "Resubstitution",
               classifier = "elasticNetGLM",
               multiViewMethod = "none",
               assayCombinations = NULL,
               nFolds = 5,
               nRepeats = 100,
               nCores = 1,
               characteristicsLabel = NULL)

{
    # Check that data is in the right format
    checkData(measurements, outcome)
    
    # Check that other variables are in the right format and fix
    nFeatures <- cleanNFeatures(nFeatures = nFeatures,
                                measurements = measurements)
    selectionMethod <- cleanSelectionMethod(selectionMethod = selectionMethod,
                                            measurements = measurements)
    classifier <- cleanClassifier(classifier = classifier,
                                  measurements = measurements)

    # Which data-types or data-views are present?
    if(is.null(characteristicsLabel)) characteristicsLabel <- "none"

    # Setup cross-validation parameters including
    crossValParams <- generateCrossValParams(nRepeats = nRepeats,
                                             nFolds = nFolds,
                                             nCores = nCores,
                                             selectionOptimisation = selectionOptimisation
    )

    # Turn text into TrainParams and TestParams objects
    modellingParams <- generateModellingParams(assayIDs = assayIDs,
                                               measurements = measurements,
                                               nFeatures = nFeatures,
                                               selectionMethod = selectionMethod,
                                               selectionOptimisation = selectionOptimisation,
                                               classifier = classifier,
                                               multiViewMethod = multiViewMethod
    )
    if(length(assayIDs) > 1 || length(assayIDs) == 1 && assayIDs != 1) assayText <- assayIDs else assayText <- NULL
    characteristics <- S4Vectors::DataFrame(characteristic = c(if(!is.null(assayText)) "Assay Name" else NULL, "Classifier Name", "Selection Name", "multiViewMethod", "characteristicsLabel"), value = c(if(!is.null(assayText)) paste(assayText, collapse = ", ") else NULL, paste(classifier, collapse = ", "),  paste(selectionMethod, collapse = ", "), multiViewMethod, characteristicsLabel))

    classifyResults <- runTests(measurements, outcome, crossValParams = crossValParams, modellingParams = modellingParams, characteristics = characteristics)
    
    fullResult <- runTest(measurements, outcome, measurements, outcome, crossValParams = crossValParams, modellingParams = modellingParams, characteristics = characteristics, .iteration = 1)

    classifyResults@finalModel <- list(fullResult$models)
    classifyResults

}






simplifyResults <- function(results, values = c("assay", "classifier", "selectionMethod", "multiViewMethod")){
    ch <- sapply(results, function(x) x@characteristics[x@characteristics$characteristic %in% values, "value"], simplify = TRUE)
    ch <- data.frame(t(ch))
    results[!duplicated(ch)]
}

#' @rdname crossValidate
#' @export
setMethod("predict", "ClassifyResult", 
          function(object, newData)
          {
              object@modellingParams@predictParams@predictor(object@finalModel[[1]], newData)
          })
