#' Automated Selection of Previously Selected Features
#' 
#' Uses the feature selection of the same cross-validation iteration of a
#' previous classification for the current classification task.
#' 
#' 
#' @aliases previousSelection previousSelection,matrix-method
#' previousSelection,DataFrame-method
#' previousSelection,MultiAssayExperiment-method
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix}, the rows are features, and the columns are samples.
#' @param classes Do not specify this variable. It is ignored and only used to
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
#'     genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
#'     genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
#'                                  c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
#'     colnames(genesMatrix) <- paste("Sample", 1:50)
#'     rownames(genesMatrix) <- paste("Gene", 1:100)                                 
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#' 
#'     CVparams <- CrossValParams(permutations = 2, folds = 2)
#'     result <- runTests(genesMatrix, classes, CVparams, ModellingParams())
#'     features(result)
#'                        
#'     # Genes 50 to 74 have differential expression in new data set.
#'     newDataset <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
#'     newDataset <- cbind(newDataset, rbind(sapply(1:25, function(sample) rnorm(49, 9, 2)),
#'                                           sapply(1:25, function(sample) rnorm(25, 14, 2)),
#'                                           sapply(1:25, function(sample) rnorm(26, 9, 2))))
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
#'     features(newerResult)
#'   #}  
#' 
#' 
#' @importFrom MultiAssayExperiment colData wideFormat
#' @importFrom S4Vectors mcols
#' @export
setGeneric("previousSelection", function(measurements, ...)
standardGeneric("previousSelection"))

setMethod("previousSelection", "matrix", 
          function(measurements, ...)
{
  previousSelection(DataFrame(t(measurements), check.names = FALSE), ...)
})

# Classes is passed around because most other selection functions need it, so it is sent from
# .doSelection but of course not used here.
setMethod("previousSelection", "DataFrame", 
          function(measurements, classes, classifyResult, minimumOverlapPercent = 80,
                   .iteration, verbose = 3)
{
  if(verbose == 3)
    message("Choosing previous features.")
  
  previousIDs <- features(classifyResult)[[.iteration]]
  if(is.character(previousIDs))
  {
    commonFeatures <- intersect(previousIDs, colnames(measurements))
    overlapPercent <- length(commonFeatures) / length(previousIDs) * 100
  } else { # A data.frame describing the data set and variable name of the chosen feature.
    keepRows <- numeric()
    varInfo <- S4Vectors::mcols(measurements) # mcols stores source information about variables.
    variable <- varInfo[, "rowname"]
    variable[is.na(variable)] <- varInfo[is.na(variable), "colname"]
    for(index in 1:length(previousIDs))
    {
      if(any(previousIDs[index, "dataset"] == varInfo[, "sourceName"] & previousIDs[index, "variable"] == variable))
        keepRows <- c(keepRows, index)
    }
    commonFeatures <- previousIDs[keepRows, ]
    overlapPercent <- nrow(commonFeatures) / nrow(previousIDs) * 100
  }
  if(overlapPercent < minimumOverlapPercent)
    signalCondition(simpleError(paste("Number of features in common between previous and current data set is lower than", minimumOverlapPercent, "percent.")))
  
  commonFeatures # Ranking isn't transferred across.
})

setMethod("previousSelection", "MultiAssayExperiment", 
          function(measurements, ...)
          {
            clinicalColumns <- colnames(MultiAssayExperiment::colData(clinicalColumns))
            dataTable <- wideFormat(measurements, colDataCols = clinicalColumns, check.names = FALSE, collapse = ':')
            S4Vectors::mcols(dataTable)[, "sourceName"] <- gsub("colDataCols", "clinical", S4Vectors::mcols(dataTable)[, "sourceName"])
            previousSelection(dataTable, ...)
          })