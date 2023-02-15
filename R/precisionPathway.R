# Basically, an ordered list of cross-validations.

setGeneric("precisionPathwayTrain", function(measurements, class, ...)
    standardGeneric("precisionPathwayTrain"))

#' @rdname precisionPathwayTrain
#' @export
setMethod("precisionPathwayTrain", "MultiAssayExperiment", 
          function(measurements, class, clinicalPredictors = NULL, ..., fixedAssays = "clinical",
                   confidenceCutoff = 0.8, minAssaySamples = 10,
                   nFeatures = 20, selectionMethod = setNames(c(NULL, rep("t-test", length(measurements))), c("clinical", names(measurements))),
                   classifier = setNames(c("elasticNetGLM", rep("randomForest", length(measurements))), c("clinical", names(measurements))),
                   nFolds = 5, nRepeats = 20, nCores = 1)
          {
              prepArgs <- list(measurements, class, clinicalPredictors)              
              extraInputs <- list(...)
              prepExtras <- numeric()
              if(length(extraInputs) > 0)
                prepExtras <- which(names(extraInputs) %in% .ClassifyRenvir[["prepareDataFormals"]])
              if(length(prepExtras) > 0)
                prepArgs <- append(prepArgs, extraInputs[prepExtras])
              measurementsAndClass <- do.call(prepareData, prepArgs)
              
             .precisionPathwayTrain(measurementsAndClass[["measurements"]], measurementsAndClass[["outcome"]],
                                    fixedAssays = fixedAssays, confidenceCutoff = confidenceCutoff,
                                    minAssaySamples = minAssaySamples, nFeatures = nFeatures,
                                    selectionMethod = selectionMethod, classifier = classifier,
                                    nFolds = nFolds, nRepeats = nRepeats, nCores = nCores)
          })

#' @rdname precisionPathwayTrain
#' @export
setMethod("precisionPathwayTrain", "list", 
          function(measurements, class, clinicalPredictors = NULL, fixedAssays = "clinical",
                   confidenceCutoff = 0.8, minAssaySamples = 10,
                   nFeatures = 20, selectionMethod = setNames(c(NULL, rep("t-test", length(measurements))), c("clinical", names(measurements))),
                   classifier = setNames(c("elasticNetGLM", rep("randomForest", length(measurements))), c("clinical", names(measurements))),
                   nFolds = 5, nRepeats = 20, nCores = 1, ...)
          {
            # One of the tables must be named "clinical".
            if (!any(names(measurements) == "clinical"))
              stop("One of the tables must be named \"clinical\".")
              
            .precisionPathwayTrain(measurements = combined_df, class = class)
          })

# Internal method which carries out all of the processing, obtaining reformatted data from the
# MultiAssayExperiment and list (of basic rectangular tables) S4 methods.
.precisionPathwayTrain <- function(measurements, class, fixedAssays = "clinical",
                   confidenceCutoff = 0.8, minAssaySamples = 10,
                   nFeatures = 20, selectionMethod = setNames(c(NULL, rep("t-test", length(measurements))), c("clinical", names(measurements))),
                   classifier = setNames(c("elasticNetGLM", rep("randomForest", length(measurements))), c("clinical", names(measurements))),
                   nFolds = 5, nRepeats = 20, nCores = 1, ...)
          {
            # Step 1: Separate the vector classes from the table of measurements, if not already separate.
            prepArgs <- list(measurements, class, clinicalPredictors)              
            extraInputs <- list(...)
            prepExtras <- numeric()
            if(length(extraInputs) > 0)
              prepExtras <- which(names(extraInputs) %in% .ClassifyRenvir[["prepareDataFormals"]])
            if(length(prepExtras) > 0)
              prepArgs <- append(prepArgs, extraInputs[prepExtras])
            
            measurementsAndClass <- do.call(prepareData, prepArgs)
            measurements <- measurementsAndClass[["measurements"]]
            class <- measurementsAndClass[["outcome"]]
            
            # Step 2: Determine all valid permutations of assays, taking into account the
            # assays to be used and which assays, if any, must be included.
            assayIDs <- unique(S4Vectors::mcols(measurements)[["assay"]])
            assaysPermutations <- .permutations(assayIDs, fixed = data.frame(seq_along(fixedAssays), fixedAssays))
            
            # Step 3: Build a classifier for each assay using all of the samples.
            modelsList <- crossValidate(measurements, class, nFeatures, selectionMethod,
                                        classifier = classifier, nFolds = nFolds,
                                        nRepeats = nRepeats, nCores = nCores)
        }