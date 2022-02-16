#' Ranking of Differential Distributions with Kullback-Leibler Distance
#' 
#' Ranks features from largest Kullback-Leibler distance between classes to
#' smallest.
#' 
#' The distance is defined as \eqn{\frac{1}{2} \times \big( \frac{(location_1 -
#' location_2)^2}{scale_1^2} + \frac{(location_1 - location_2)^2}{scale_2^2} +
#' \frac{scale_2^2}{scale_1^2} + \frac{scale_1^2}{scale_2^2}\big)}{0.5 *
#' ((location1 - location2)^2 / scale1^2 + (location1 - location2)^2 / scale2^2
#' + scale1^2 / scale2^2 + scale2^2 / scale1^2)}
#' 
#' The subscripts denote the group which the parameter is calculated for.
#' 
#' @aliases KullbackLeiblerRanking KullbackLeiblerRanking,matrix-method
#' KullbackLeiblerRanking,DataFrame-method
#' KullbackLeiblerRanking,MultiAssayExperiment-method
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
#' \code{\link{getLocationsAndScales}}.
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
#'   genesMatrix <- sapply(1:20, function(feature)
#'                               {
#'                                 randomMeans <- sample(c(8, 12), 25, replace = TRUE)
#'                                 c(rnorm(25, randomMeans, 1), rnorm(25, 10, 1))
#'                               }
#'                        )
#'   genesMatrix <- cbind(genesMatrix, sapply(1:80, function(feature) rnorm(50, 10, 1)))
#'   classes <- factor(rep(c("Poor", "Good"), each = 25))
#'   
#'   ranked <- KullbackLeiblerRanking(genesMatrix, classes)
#'   head(ranked)
#' 
#' @export
setGeneric("KullbackLeiblerRanking", function(measurementsTrain, ...)
           standardGeneric("KullbackLeiblerRanking"))

# Matrix of numeric measurements.
setMethod("KullbackLeiblerRanking", "matrix", function(measurementsTrain, classesTrain, ...)
{
  KullbackLeiblerRanking(DataFrame(measurementsTrain, check.names = FALSE), classesTrain, ...)
})

setMethod("KullbackLeiblerRanking", "DataFrame", # Sample information data or one of the other inputs, transformed.
          function(measurementsTrain, classesTrain, ..., verbose = 3)
{
  splitDataset <- .splitDataAndOutcomes(measurementsTrain, classesTrain)
  measurementsTrain <- splitDataset[["measurements"]]
  
  if(verbose == 3)
    message("Selecting features by Kullback-Leibler divergence.")

  oneClassMeasurements <- measurementsTrain[classesTrain == levels(classesTrain)[1], ]
  otherClassMeasurements <- measurementsTrain[classesTrain == levels(classesTrain)[2], ]
  oneClassDistribution <- getLocationsAndScales(oneClassMeasurements, ...)
  otherClassDistribution <- getLocationsAndScales(otherClassMeasurements, ...)
  locationDifference <- oneClassDistribution[[1]] - otherClassDistribution[[1]]
  divergence <- 1/2 * (locationDifference^2 / ((oneClassDistribution[[2]])^2) +
                         locationDifference^2 / ((otherClassDistribution[[2]])^2) +
                         ((oneClassDistribution[[2]])^2) / ((otherClassDistribution[[2]])^2) +
                         ((otherClassDistribution[[2]])^2) / ((oneClassDistribution[[2]])^2))

  if(!is.null(S4Vectors::mcols(measurementsTrain)))
    S4Vectors::mcols(measurementsTrain)[order(divergence, decreasing = TRUE), ]
  else
    colnames(measurementsTrain)[order(divergence, decreasing = TRUE)]
  
})

# One or more omics data sets, possibly with sample information data.
setMethod("KullbackLeiblerRanking", "MultiAssayExperiment",
          function(measurementsTrain, targets = names(measurementsTrain), classesTrain, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurementsTrain, targets, classesTrain)
  measurementsTrain <- tablesAndClasses[["dataTable"]]
  classesTrain <- tablesAndClasses[["outcomes"]]

  if(ncol(dataTable) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    KullbackLeiblerRanking(measurementsTrain, classesTrain, ...)
})