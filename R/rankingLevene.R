#' Selection of Differential Variability with Levene Statistic
#' 
#' Ranks features by largest Levene statistic.
#' 
#' Levene's statistic for unequal variance between groups is a robust version
#' of Bartlett's statistic.
#' 
#' @aliases leveneRanking leveneRanking,matrix-method
#' leveneRanking,DataFrame-method leveneRanking,MultiAssayExperiment-method
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
#' If of type \code{\link{DataFrame}} or \code{\link{MultiAssayExperiment}}, the data set is subset
#' to only those features of type \code{numeric}.
#' @param classesTrain Either a vector of class labels of class \code{\link{factor}}
#' of the same length as the number of samples in \code{measurementsTrain} or if the
#' measurements are of class \code{DataFrame} a character vector of length 1
#' containing the column name in \code{measurement} is also permitted.
#' @param targets If \code{measurementsTrain} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"sampleInfo"} is also a valid value
#' and specifies that numeric variables from the sample information table will be
#' used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return A vector or data frame (if \code{MultiAssayExperiment} input) of
#' features, from the most promising features in the first position to the
#' least promising feature in the last position.
#' @author Dario Strbenac
#' @examples
#' 
#'   # First 20 features have bimodal distribution for Poor class.
#'   # Other 80 features have normal distribution for both classes.
#'   set.seed(1984)
#'   genesMatrix <- sapply(1:20, function(feature)
#'                               {
#'                                 randomMeans <- sample(c(8, 12), 25, replace = TRUE)
#'                                 c(rnorm(25, randomMeans, 1), rnorm(25, 10, 1))
#'                               }
#'                        )
#'   genesMatrix <- cbind(genesMatrix, sapply(1:80, function(feature) rnorm(50, 10, 1)))
#'   classes <- factor(rep(c("Poor", "Good"), each = 25))
#'   
#'   ranked <- leveneRanking(genesMatrix, classes)
#'   head(ranked)
#' 
#' @usage NULL
#' @export
setGeneric("leveneRanking", function(measurementsTrain, ...)
           standardGeneric("leveneRanking"))

# Matrix of numeric measurements.
#' @rdname leveneRanking
#' @export
setMethod("leveneRanking", "matrix", function(measurementsTrain, classesTrain, ...)
{
  leveneRanking(S4Vectors::DataFrame(measurementsTrain, check.names = FALSE), classesTrain, ...)
})

#' @rdname leveneRanking
#' @export
setMethod("leveneRanking", "DataFrame", # Sample information data or one of the other inputs, transformed.
          function(measurementsTrain, classesTrain, verbose = 3)
{
  splitDataset <- .splitDataAndOutcomes(measurementsTrain, classesTrain)
  measurementsTrain <- splitDataset[["measurements"]]
  
  if(!requireNamespace("car", quietly = TRUE))
    stop("The package 'car' could not be found. Please install it.")            
  if(verbose == 3)
    message("Calculating Levene statistic.")

  pValues <- apply(measurementsTrain, 2, function(featureColumn)
             car::leveneTest(featureColumn, classesTrain)[["Pr(>F)"]][1])
  
  order(pValues) # From smallest to largest.
})

# One or more omics data sets, possibly with sample information data.
#' @rdname leveneRanking
#' @export
setMethod("leveneRanking", "MultiAssayExperiment",
          function(measurementsTrain, targets = names(measurementsTrain), classesTrain, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurementsTrain, targets, classesTrain)
  measurementsTrain <- tablesAndClasses[["dataTable"]]
  classesTrain <- tablesAndClasses[["outcomes"]]
  
  leveneRanking(measurementsTrain, classesTrain, ...)
})