#' Ranking of Differential Distributions with Likelihood Ratio Statistic
#' 
#' Ranks features from largest difference of log likelihoods (null hypothesis -
#' alternate hypothesis) to smallest.
#' 
#' Likelihood ratio test of null hypothesis that the location and scale are the
#' same for both groups, and an alternate hypothesis that is specified by
#' parameters. The location and scale of features is calculated by
#' \code{\link{getLocationsAndScales}}. The distribution fitted to the data is
#' the normal distribution.
#' 
#' @aliases likelihoodRatioRanking likelihoodRatioRanking,matrix-method
#' likelihoodRatioRanking,DataFrame-method
#' likelihoodRatioRanking,MultiAssayExperiment-method
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
#' If of type \code{\link{DataFrame}} or \code{\link{MultiAssayExperiment}}, the data set is subset
#' to only those features of type \code{numeric}.
#' @param classesTrain Either a vector of class labels of class \code{\link{factor}}
#' of the same length as the number of samples in \code{measurementsTrain} or if the
#' measurements are of class \code{DataFrame} a character vector of length 1
#' containing the column name in \code{measurement} is also permitted. Not used
#' if \code{measurementsTrain} is a \code{MultiAssayExperiment} object.
#' @param targets If \code{measurementsTrain} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"sampleInfo"} is also a valid value
#' and specifies that numeric variables from the sample information data table will be
#' used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method or options which are accepted by the function
#' \code{\link{getLocationsAndScales}}.
#' @param alternative Default: \code{c("different", "different")}. A vector of
#' length 2.  The first element specifies the location of the alternate
#' hypothesis.  The second element specifies the scale of the alternate
#' hypothesis.  Valid values in each element are "same" or "different".
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
#' 
#'   genesMatrix <- sapply(1:20, function(feature)
#'                               {
#'                                 randomMeans <- sample(c(8, 12), 25, replace = TRUE)
#'                                 c(rnorm(25, randomMeans, 1), rnorm(25, 10, 1))
#'                               }
#'                        )
#'   genesMatrix <- cbind(genesMatrix, sapply(1:80, function(feature) rnorm(50, 10, 1)))
#'   classes <- factor(rep(c("Poor", "Good"), each = 25))
#' 
#'   ranked <- likelihoodRatioRanking(genesMatrix, classes)
#'   head(ranked)
#'
#' @usage NULL
#' @export
setGeneric("likelihoodRatioRanking", function(measurementsTrain, ...)
           standardGeneric("likelihoodRatioRanking"))

# Matrix of numeric measurements.
#' @rdname likelihoodRatioRanking
#' @export
setMethod("likelihoodRatioRanking", "matrix", function(measurementsTrain, classesTrain, ...)
{
  likelihoodRatioRanking(S4Vectors::DataFrame(measurementsTrain, check.names = FALSE), classesTrain, ...)
})

#' @rdname likelihoodRatioRanking
#' @export
setMethod("likelihoodRatioRanking", "DataFrame", # Sample information data or one of the other inputs, transformed.
          function(measurementsTrain, classesTrain, alternative = c(location = "different", scale = "different"),
                   ..., verbose = 3)
{
  splitDataset <- .splitDataAndOutcomes(measurementsTrain, classesTrain)
  measurementsTrain <- splitDataset[["measurements"]]

  if(verbose == 3)
    message("Ranking features by likelihood ratio test statistic.")

  allDistribution <- getLocationsAndScales(measurementsTrain, ...)
  logLikelihoodRatios <- unlist(mapply(function(featureMeasurements, scale, location)
  sum(dnorm(featureMeasurements, scale, location, log = TRUE)),
  measurementsTrain, allDistribution[[1]], allDistribution[[2]])) -
  rowSums(sapply(levels(classesTrain), function(class)
  {
    classMeasurements <- measurementsTrain[which(classesTrain == class), ]
    classDistribution <- getLocationsAndScales(classMeasurements, ...)
    
    unlist(mapply(function(featureMeasurements, scale, location)
    sum(dnorm(featureMeasurements, scale, location, log = TRUE)),
    classMeasurements,
    switch(alternative[["location"]], same = allDistribution[[1]], different = classDistribution[[1]]),
    switch(alternative[["scale"]], same = allDistribution[[2]], different = classDistribution[[2]])))    
  }))
  
  if(!is.null(S4Vectors::mcols(measurementsTrain)))
    S4Vectors::mcols(measurementsTrain)[order(logLikelihoodRatios), ]
  else
    colnames(measurementsTrain)[order(logLikelihoodRatios)]
})

# One or more omics data sets, possibly with sample information data.
#' @rdname likelihoodRatioRanking
#' @export
setMethod("likelihoodRatioRanking", "MultiAssayExperiment",
          function(measurementsTrain, targets = names(measurementsTrain), classesTrain, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurementsTrain, targets, classesTrain)
  measurementsTrain <- tablesAndClasses[["dataTable"]]
  classesTrain <- tablesAndClasses[["outcomes"]]

  if(ncol(measurementsTrain) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
  likelihoodRatioRanking(measurementsTrain, classesTrain, ...)
})