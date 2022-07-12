#' Ranking of Pairs of Features that are Different Between Classes
#' 
#' Ranks pre-specified pairs of features by the largest difference of the sum
#' of measurement differences over all samples within a class.
#' 
#' Instead of considering whether one feature in a pair of features is
#' consistently lower or higher than the other in the pair, this method takes
#' the sum of differences across all samples within a class, to prevent ties in
#' the ranking of pairs of features.
#' 
#' @aliases pairsDifferencesRanking pairsDifferencesRanking,matrix-method
#' pairsDifferencesRanking,DataFrame-method
#' pairsDifferencesRanking,MultiAssayExperiment-method
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
#' If of type \code{\link{DataFrame}} or \code{\link{MultiAssayExperiment}}, the data set is subset
#' to only those features of type \code{numeric}.
#' @param classesTrain Either a vector of class labels of class \code{\link{factor}}
#' of the same length as the number of samples in \code{measurementsTrain} or if the
#' measurements are of class \code{DataFrame} a character vector of length 1
#' containing the column name in \code{measurement} is also permitted.
#' @param featurePairs An S4 object of type \code{\link{Pairs}} containing
#' feature identifiers to calculate the sum of differences within each class
#' for.
#' @param target If \code{measurementsTrain} is a \code{MultiAssayExperiment}, the
#' name of the data table to be used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return A \code{\link{Pairs}} object, from the most promising feature pair
#' in the first position to the least promising feature pair in the last
#' position.
#' @author Dario Strbenac
#' @seealso \code{\link{kTSPclassifier}} for a classifier which makes use of
#' the pairs of selected features in classification.
#' @references Simple decision rules for classifying human cancers from gene
#' expression profiles, Aik C Tan, Daniel Q Naiman, Lei Xu, Raimond L. Winslow
#' and Donald Geman, 2005, \emph{Bioinformatics}, Volume 21 Issue 20,
#' \url{https://academic.oup.com/bioinformatics/article/21/20/3896/203010}.
#' @examples
#' 
#'   featurePairs <- Pairs(c('A', 'C'), c('B', 'C'))
#'                            
#'   # Difference in differences for features A and C between classes.                                           
#'   measurements <- matrix(c(9.9, 10.5, 10.1, 10.9, 11.0, 6.6, 7.7, 7.0, 8.1, 6.5,
#'                            8.5, 10.5, 12.5, 10.5, 9.5, 8.5, 10.5, 12.5, 10.5, 9.5,
#'                            6.6, 7.7, 7.0, 8.1, 6.5, 11.2, 11.0, 11.1, 11.4, 12.0,
#'                            8.1, 10.6, 7.4, 7.1, 10.4, 6.1, 7.3, 2.7, 11.0, 9.1,
#'                            round(rnorm(60, 8, 1), 1)), nrow = 10)
#'   classes <- factor(rep(c("Good", "Poor"), each = 5))
#'                          
#'   rownames(measurements) <- paste("Patient", 1:10)
#'   colnames(measurements) <- LETTERS[1:10]
#'   
#'   pairsDifferencesRanking(measurements, classes, featurePairs = featurePairs)
#' 
#' @usage NULL
#' @export
setGeneric("pairsDifferencesRanking", function(measurementsTrain, ...)
           standardGeneric("pairsDifferencesRanking"))

#' @rdname pairsDifferencesRanking
#' @export
setMethod("pairsDifferencesRanking", "matrix", # Matrix of numeric measurements.
          function(measurementsTrain, classesTrain, featurePairs = NULL, ...)
{
  if(is.null(featurePairs))
    stop("No feature pairs provided but some must be.")
  if(!"Pairs" %in% class(featurePairs))
    stop("'featurePairs' must be of type Pairs.")
            
  pairsDifferencesRanking(S4Vectors::DataFrame(measurementsTrain, check.names = FALSE), classesTrain, featurePairs, ...)
})

#' @rdname pairsDifferencesRanking
#' @export
setMethod("pairsDifferencesRanking", "DataFrame",
          function(measurementsTrain, classesTrain, featurePairs = NULL, verbose = 3)
{
  if(is.null(featurePairs))
    stop("No feature pairs provided but some must be.")
  if(!"Pairs" %in% class(featurePairs))
    stop("'featurePairs' must be of type Pairs.")
            
  splitDataset <- .splitDataAndOutcomes(measurementsTrain, classesTrain)
  measurementsTrain <- splitDataset[["measurements"]]
  
  suppliedPairs <- length(featurePairs)
  keepPairs <- S4Vectors::first(featurePairs) %in% colnames(measurementsTrain) & S4Vectors::second(featurePairs) %in% colnames(measurementsTrain)
  featurePairs <- featurePairs[keepPairs]
  
  if(verbose == 3)
    message(suppliedPairs, " pairs input and ", length(featurePairs), " pairs remain after filtering based on data set row names.")
  
  if(verbose == 3)
    message("Selecting pairs of features with consistent differences.")

  oneClassTraining <- which(classesTrain == levels(classesTrain)[1])
  otherClassTraining <- which(classesTrain == levels(classesTrain)[2])
  oneClassMeasurements <- measurementsTrain[oneClassTraining, ]
  otherClassMeasurements <- measurementsTrain[otherClassTraining, ]

  numerator <- as.matrix(oneClassMeasurements[, S4Vectors::first(featurePairs)])
  denominator <- as.matrix(oneClassMeasurements[, S4Vectors::second(featurePairs)])
  oneClassDifferences <- colMeans(numerator - denominator)
  
  numerator <- as.matrix(otherClassMeasurements[, S4Vectors::first(featurePairs)])
  denominator <- as.matrix(otherClassMeasurements[, S4Vectors::second(featurePairs)])
  otherClassDifferences <- colMeans(numerator - denominator)
  
  pairsClassDifferences <- otherClassDifferences - oneClassDifferences
  
  order(abs(pairsClassDifferences), decreasing = TRUE)
})

# One or more omics data sets, possibly with sample information data.
#' @rdname pairsDifferencesRanking
#' @export
setMethod("pairsDifferencesRanking", "MultiAssayExperiment",
          function(measurementsTrain, target = names(measurementsTrain)[1], classesTrain, featurePairs = NULL, ...)
{
  if(is.null(featurePairs))
    stop("No feature pairs provided but some must be.")
  if(!"Pairs" %in% class(featurePairs))
    stop("'featurePairs' must be of type Pairs.")         
            
  tablesAndClasses <- .MAEtoWideTable(measurementsTrain, target, classesTrain)
  measurementsTrain <- tablesAndClasses[["dataTable"]]
  classesTrain <- tablesAndClasses[["outcomes"]]            
  pairsDifferencesRanking(measurementsTrain, classesTrain, featurePairs, ...)
})