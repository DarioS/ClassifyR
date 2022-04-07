#' Automated Selection of Previously Selected Features
#' 
#' Uses the feature selection of the same cross-validation iteration of a
#' previous classification for the current classification task.
#' 
#' 
#' @aliases previousSelection previousSelection,matrix-method
#' previousSelection,DataFrame-method
#' previousSelection,MultiAssayExperiment-method
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
#' \code{matrix}, the rows are features, and the columns are samples.
#' @param classesTrain Do not specify this variable. It is ignored and only used to
#' create consistency of formal parameters with other feature selection
#' methods.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method.
#' @param classifyResult An existing classification result from which to take
#' the feature selections from.
#' @param minimumOverlapPercent If at least this many selected features can't
#' be identified in the current data set, then the selection stops with an
#' error.
#' @param .iteration Do not specify this variable. It is set by
#' \code{\link{runTests}} if this function is being repeatedly called by
#' \code{runTests}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return A vector of feature indices, from the most promising features in the
#' first position to the least promising feature in the last position.
#' @author Dario Strbenac
#' @examples
#' 
#'   #if(require(sparsediscrim))
#'   #{
#'     # Genes 76 to 100 have differential expression.
#'     genesMatrix <- sapply(1:100, function(sample) rnorm(25, 9, 0.3))
#'     genesMatrix <- rbind(genesMatrix, t(sapply(1:25, function(sample)
#'                                       c(rnorm(75, 9, 0.3), rnorm(25, 14, 0.3)))))
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     rownames(genesMatrix) <- paste("Sample", 1:50)
#'     colnames(genesMatrix) <- paste("Gene", 1:100)                          
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#' 
#'     CVparams <- CrossValParams(permutations = 2, folds = 2)
#'     result <- runTests(genesMatrix, classes, CVparams, ModellingParams())
#'     chosenFeatureNames(result)
#'                        
#'     # Genes 50 to 74 have differential expression in new data set.
#'     
#'     newDataset <- sapply(1:100, function(sample) c(rnorm(25, 9, 0.3)))
#'     newDataset <- rbind(newDataset, cbind(sapply(1:49, function(sample) rnorm(25, 9, 0.3)),
#'                                           sapply(1:25, function(sample) rnorm(25, 14, 0.3)),
#'                                           sapply(1:26, function(sample) rnorm(25, 9, 0.3))))
#'                                           
#'     rownames(newDataset) <- rownames(genesMatrix)
#'     colnames(newDataset) <- colnames(genesMatrix)
#'     
#'     selPars <- SelectParams(previousSelection, intermediate = ".iteration", classifyResult = result)
#'     previousParams <- ModellingParams(selectParams = selPars)                                       
#'     newerResult <- runTests(newDataset, classes, CVparams, previousParams)
#'                                      
#'     # However, only genes 76 to 100 are chosen, because the feature selections are
#'     # carried over from the first cross-validated classification.
#'     chosenFeatureNames(newerResult)
#'   #}  
#' 
#' 
#' @importFrom MultiAssayExperiment colData wideFormat
#' @importFrom S4Vectors mcols
#' @usage NULL
#' @export
setGeneric("previousSelection", function(measurementsTrain, ...)
standardGeneric("previousSelection"))

#' @rdname previousSelection
#' @export
setMethod("previousSelection", "matrix", 
          function(measurementsTrain, ...)
{
  previousSelection(S4Vectors::DataFrame(measurementsTrain, check.names = FALSE), ...)
})

# Classes is passed around because most other selection functions need it, so it is sent from
# .doSelection but of course not used here.
#' @rdname previousSelection
#' @export
setMethod("previousSelection", "DataFrame", 
          function(measurementsTrain, classesTrain, classifyResult, minimumOverlapPercent = 80,
                   .iteration, verbose = 3)
{
  if(verbose == 3)
    message("Choosing previous features.")

  previousIDs <- chosenFeatureNames(classifyResult)[[.iteration]]
  if(is.character(previousIDs))
  {
    commonFeatures <- intersect(previousIDs, colnames(measurementsTrain))
    overlapPercent <- length(commonFeatures) / length(previousIDs) * 100
  } else { # A data.frame describing the data set and variable name of the chosen feature.
    featuresIDs <- do.call(paste, S4Vectors::mcols(measurementsTrain)[, c("dataset", "feature")])
    selectedIDs <-  do.call(paste, previousIDs)
    selectedColumns <- match(selectedIDs, featuresIDs)
    commonFeatures <- sum(!is.na(selectedColumns))
    overlapPercent <- commonFeatures / nrow(previousIDs) * 100
  }
  if(overlapPercent < minimumOverlapPercent)
    signalCondition(simpleError(paste("Number of features in common between previous and current data set is lower than", minimumOverlapPercent, "percent.")))
  
  S4Vectors::mcols(measurementsTrain)[selectedColumns, ] # Each row is about one column.
})

#' @rdname previousSelection
#' @export
setMethod("previousSelection", "MultiAssayExperiment", 
          function(measurementsTrain, ...)
          {
            sampleInfoColumns <- colnames(MultiAssayExperiment::colData(sampleInfoColumns))
            dataTable <- wideFormat(measurementsTrain, colDataCols = sampleInfoColumns, check.names = FALSE, collapse = ':')
            S4Vectors::mcols(dataTable)[, "sourceName"] <- gsub("colDataCols", "sampleInfo", S4Vectors::mcols(dataTable)[, "sourceName"])
            previousSelection(dataTable, ...)
          })