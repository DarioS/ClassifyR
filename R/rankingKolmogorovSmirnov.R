#' Ranking of Differential Distributions with Kolmogorov-Smirnov Distance
#' 
#' Ranks features from largest Kolmogorov-Smirnov distance to smallest.
#' 
#' @aliases KolmogorovSmirnovRanking KolmogorovSmirnovRanking,matrix-method
#' KolmogorovSmirnovRanking,DataFrame-method
#' KolmogorovSmirnovRanking,MultiAssayExperiment-method
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
#' and specifies that numeric variables from the sample information data table will be
#' used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method or options which are accepted by the function
#' \code{\link{ks.test}}.
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
#'   genesMatrix <- sapply(1:25, function(sample)
#'                               {
#'                                 randomMeans <- sample(c(8, 12), 20, replace = TRUE)
#'                                 c(rnorm(20, randomMeans, 1), rnorm(80, 10, 1))
#'                               }
#'                        )
#'   genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample) rnorm(100, 10, 1)))
#'   rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
#'   classes <- factor(rep(c("Poor", "Good"), each = 25))
#' 
#'   ranked <- KolmogorovSmirnovRanking(genesMatrix, classes)
#'   head(ranked)                           
#' 
#' @export
setGeneric("KolmogorovSmirnovRanking", function(measurementsTrain, ...)
           standardGeneric("KolmogorovSmirnovRanking"))

# Matrix of numeric measurements.
setMethod("KolmogorovSmirnovRanking", "matrix", function(measurementsTrain, classesTrain, ...)
{
  KolmogorovSmirnovRanking(DataFrame(measurementsTrain, check.names = FALSE), classesTrain, ...)
})

setMethod("KolmogorovSmirnovRanking", "DataFrame", # Sample information data or one of the other inputs, transformed.
          function(measurementsTrain, classesTrain, ..., verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurementsTrain, classesTrain)
  measurementsTrain <- splitDataset[["measurements"]]
  
  if(verbose == 3)
    message("Ranking features by Kolmogorov Smirnov distance between classes.")

  oneClass <- classesTrain == levels(classesTrain)[1]
  otherClass <- classesTrain == levels(classesTrain)[2]
  KSdistance <- apply(measurementsTrain, 2, function(featureColumn)
                      stats::ks.test(featureColumn[oneClass], featureColumn[otherClass], ...)[["statistic"]])

  if(!is.null(S4Vectors::mcols(measurementsTrain)))
    S4Vectors::mcols(measurementsTrain)[order(KSdistance, decreasing = TRUE), ]
  else
    colnames(measurementsTrain)[order(KSdistance, decreasing = TRUE)]
})

# One or more omics data sets, possibly with sample information data.
setMethod("KolmogorovSmirnovRanking", "MultiAssayExperiment",
function(measurementsTrain, targets = names(measurementsTrain), classesTrain, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurementsTrain, targets, classesTrain)
  measurementsTrain <- tablesAndClasses[["dataTable"]]
  classesTrain <- tablesAndClasses[["outcomes"]]
            
  if(ncol(dataTable) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    KolmogorovSmirnovRanking(measurementsTrain, classesTrain, ...)
})