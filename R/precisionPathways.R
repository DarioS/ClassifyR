#' Precision Pathways for Sample Prediction Based on Prediction Confidence.
#' 
#' Precision pathways allows the evaluation of various permutations of multiomics or multiview data.
#' Samples are predicted by a particular assay if they were consistently predicted as a particular class
#' during cross-validation. Otherwise, they are passed onto subsequent assays/tiers for prediction. Balanced accuracy
#' is used to evaluate overall prediction performance and sample-specific accuracy for individual-level evaluation.
#'
#' @param measurements Either a \code{\link{MultiAssayExperiment}} or a list of the basic tabular objects containing the data.
#' @param class Same as \code{measurements} but only training samples. IF \code{measurements} is a \code{list}, may also be
#' a vector of classes.
#' @param clinicalPredictors Default: \code{NULL}. Must be a character vector of clinical features to use in modelling. This allows avoidance of things like sample IDs,
#' sample acquisition dates, etc. which are not relevant for outcome prediction.
#' @param maxMissingProp Default: 0.0. A proportion less than 1 which is the maximum
#' tolerated proportion of missingness for a feature to be retained for modelling.
#' @param topNvariance Default: NULL. An integer number of most variable features per assay to subset to.
#' Assays with less features won't be reduced in size.
#' @param fixedAssays A character vector of assay names specifying any assays which must be at the
#' beginning of the pathway.
#' @param confidenceCutoff The minimum confidence of predictions for a sample to be predicted by a particular issue
#' . If a sample was predicted to belong to a particular class a proportion \eqn{p} times, then the confidence is \eqn{2 \times |p - 0.5|}.
#' @param minAssaySamples An integer specifying the minimum number of samples a tier may have. If a subsequent tier
#' would have less than this number of samples, the samples are incorporated into the current tier.
#' @param nFeatures Default: 20. The number of features to consider during feature selection, if feature selection is done.
#' @param selectionMethod A named character vector of feature selection methods to use for the assays, one for each. The names must correspond to names of \code{measurements}.
#' @param classifier A named character vector of modelling methods to use for the assays, one for each. The names must correspond to names of \code{measurements}.
#' @param nFolds A numeric specifying the number of folds to use for cross-validation.
#' @param nRepeats A numeric specifying the the number of repeats or permutations to use for cross-validation.
#' @param nCores A numeric specifying the number of cores used if the user wants to use parallelisation.
#' @param pathways A set of pathways created by \code{precisionPathwaysTrain} to be used for predicting on a new data set.
#' @rdname precisionPathways
#' @return An object of class \code{PrecisionPathways} which is basically a named list that other plotting and
#' tabulating functions can use.
#' @examples
#' # To be determined.

#' @usage NULL
setGeneric("precisionPathwaysTrain", function(measurements, class, ...)
    standardGeneric("precisionPathwaysTrain"))

#' @rdname precisionPathways
#' @export
setMethod("precisionPathwaysTrain", "MultiAssayExperimentOrList", 
          function(measurements, class, clinicalPredictors = NULL, maxMissingProp = 0.0, topNvariance = NULL,
                   fixedAssays = "clinical", confidenceCutoff = 0.8, minAssaySamples = 10,
                   nFeatures = 20, selectionMethod = setNames(c("none", rep("t-test", length(measurements))), c("clinical", names(measurements))),
                   classifier = setNames(c("elasticNetGLM", rep("randomForest", length(measurements))), c("clinical", names(measurements))),
                   nFolds = 5, nRepeats = 20, nCores = 1)
          {
            if(is.list(measurements)) # Ensure plain list has clinical data.
            {
              # One of the tables must be named "clinical".
              if (!any(names(measurements) == "clinical"))
                stop("One of the tables must be named \"clinical\".")
            }
            prepArgs <- list(measurements, outcomeColumns = class, clinicalPredictors = clinicalPredictors,
                             maxMissingProp = maxMissingProp, topNvariance = topNvariance)
            measurementsAndClass <- do.call(prepareData, prepArgs)
              
            .precisionPathwaysTrain(measurementsAndClass[["measurements"]], measurementsAndClass[["outcome"]],
                                   fixedAssays = fixedAssays, confidenceCutoff = confidenceCutoff,
                                   minAssaySamples = minAssaySamples, nFeatures = nFeatures,
                                   selectionMethod = selectionMethod, classifier = classifier,
                                   nFolds = nFolds, nRepeats = nRepeats, nCores = nCores)
          })

# Internal method which carries out all of the processing, obtaining reformatted data from the
# MultiAssayExperiment and list (of basic rectangular tables) S4 methods.
.precisionPathwaysTrain <- function(measurements, class, fixedAssays = "clinical",
                   confidenceCutoff = 0.8, minAssaySamples = 10,
                   nFeatures = 20, selectionMethod = setNames(c(NULL, rep("t-test", length(measurements))), c("clinical", names(measurements))),
                   classifier = setNames(c("elasticNetGLM", rep("randomForest", length(measurements))), c("clinical", names(measurements))),
                   nFolds = 5, nRepeats = 20, nCores = 1)
          {
            # Step 1: Determine all valid permutations of assays, taking into account the
            # assays to be used and which assays, if any, must be included.
            assayIDs <- unique(S4Vectors::mcols(measurements)[["assay"]])
            assaysPermutations <- .permutations(assayIDs, fixed = data.frame(seq_along(fixedAssays), fixedAssays))
            permutationIDs <- apply(assaysPermutations, 2, function(permutation) paste(permutation, collapse = '-'))
            
            # Step 2: Build a classifier for each assay using all of the samples.
            modelsList <- crossValidate(measurements, class, nFeatures, selectionMethod,
                                        classifier = classifier, nFolds = nFolds,
                                        nRepeats = nRepeats, nCores = nCores)
            modelsList <- lapply(modelsList, calcCVperformance, "Sample Accuracy") # Add sample accuracy, which can be subset later.

            # Step 3: Loop over each pathway and each assay in order to determine which samples are used at that level
            # and which are passed onwards.
            precisionPathways <- lapply(as.data.frame(assaysPermutations), function(permutation)
            {
              assaysProcessed <- character()
              samplesUsed <- character()
              individualsTableAll <- S4Vectors::DataFrame()
              tierTableAll <- S4Vectors::DataFrame()
              breakEarly = FALSE
              for(assay in permutation)
              {
                # Step 3a: Identify all samples which are consistently predicted.
                modelIndex <- match(assay, assayIDs)
                allPredictions <- predictions(modelsList[[modelIndex]])
                allSampleIDs <- sampleNames(modelsList[[modelIndex]])
                predictionsSamplesCounts <- table(allPredictions[, "sample"], allPredictions[, "class"])
                confidences <- 2 * abs(predictionsSamplesCounts[, 1] / rowSums(predictionsSamplesCounts) - 0.5)
                sampleIDsUse <- names(confidences)[confidences > confidenceCutoff]
                
                # Check if too few samples left for next round. Include them in this round, if so.
                remainingIDs <- setdiff(allSampleIDs, c(samplesUsed, sampleIDsUse))
                if(length(remainingIDs) < minAssaySamples)
                {
                  sampleIDsUse <- c(sampleIDsUse, remainingIDs)
                  breakEarly = TRUE
                }
                
                predictionsSamplesCounts <- predictionsSamplesCounts[sampleIDsUse, ]
                
                # Step 3b: Individuals predictions and sample-wise accuracy, tier-wise error.
                maxVotes <- apply(predictionsSamplesCounts, 1, function(sample) which.max(sample))
                predictedClasses <- factor(colnames(predictionsSamplesCounts)[maxVotes],
                                           levels = colnames(predictionsSamplesCounts))    
                individualsTable <- S4Vectors::DataFrame(Tier = assay,
                                                         `Sample ID` = sampleIDsUse,
                                                         `Predicted` = predictedClasses,
                                                         `Accuracy` = performance(modelsList[[modelIndex]])[["Sample Accuracy"]][sampleIDsUse],
                                                         check.names = FALSE)
                knownClasses <- actualOutcome(modelsList[[modelIndex]])[match(sampleIDsUse, allSampleIDs)]
                balancedAccuracy <- calcExternalPerformance(knownClasses, predictedClasses)
                tierTable <- S4Vectors::DataFrame(Tier = assay,
                                                  `Balanced Accuracy` = balancedAccuracy, check.names = FALSE)
                
                assaysProcessed <- c(assaysProcessed, assay)
                individualsTableAll <- rbind(individualsTableAll, individualsTable)
                tierTableAll <- rbind(tierTableAll, tierTable)
                samplesUsed <- c(samplesUsed, sampleIDsUse)
                
                if(breakEarly == TRUE) break
              }
              pathwayString <- paste(assaysProcessed, collapse = '-')
              parameters = list(confidenceCutoff = confidenceCutoff, minAssaySamples = minAssaySamples)
              list(models = modelsList, parameters = parameters, pathway = pathwayString,
                  individuals = individualsTableAll, tiers = tierTableAll)
            })
            
            class(precisionPathways) <- "PrecisionPathways"
            names(precisionPathways) <- sapply(precisionPathways, "[[", "pathway")
            precisionPathways
}

# A nice print method to avoid flooding the screen with lots of tables
# when result is shown in console.
print.PrecisionPathways <- function(x)
{
  cat("An object of class 'PrecisionPathways'.\n")
  cat("Pathways:\n")
  cat(paste(names(x), collapse = '\n'))
}

#' @usage NULL
setGeneric("precisionPathwaysPredict", function(pathways, measurements, class, ...)
    standardGeneric("precisionPathwaysPredict"))

#' @rdname precisionPathways
#' @export
setMethod("precisionPathwaysPredict", "MultiAssayExperimentOrList", 
          function(pathways, measurements, class)
          {
            if(is.list(measurements)) # Ensure plain list has clinical data.
            {
              # One of the tables must be named "clinical".
              if (!any(names(measurements) == "clinical"))
                stop("One of the tables must be named \"clinical\".")
            }
            prepArgs <- list(measurements, outcomeColumns = class)
            measurementsAndClass <- do.call(prepareData, prepArgs)
              
            .precisionPathwaysPredict(pathways, measurementsAndClass[["measurements"]], measurementsAndClass[["outcome"]])
          })

.precisionPathwaysPredict <- function(pathways, measurements, class)
{
  # To do.
}