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
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix}, the rows are features, and the columns are samples.
#' @param classes Either a vector of class labels of class \code{\link{factor}}
#' of the same length as the number of samples in \code{measurements} or if the
#' measurements are of class \code{DataFrame} a character vector of length 1
#' containing the column name in \code{measurement} is also permitted. Not used
#' if \code{measurements} is a \code{MultiAssayExperiment} object.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"clinical"} is also a valid value
#' and specifies that numeric variables from the clinical data table will be
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
#'   ranked <- likelihoodRatioRanking(genesMatrix, classes)
#'   head(ranked)
#' 
#' @export
setGeneric("likelihoodRatioRanking", function(measurements, ...)
           standardGeneric("likelihoodRatioRanking"))

# Matrix of numeric measurements.
setMethod("likelihoodRatioRanking", "matrix", function(measurements, classes, ...)
{
  likelihoodRatioRanking(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("likelihoodRatioRanking", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, alternative = c(location = "different", scale = "different"),
                   ..., verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")

  if(verbose == 3)
    message("Ranking features by likelihood ratio test statistic.")

  allDistribution <- getLocationsAndScales(measurements, ...)
  logLikelihoodRatios <- unlist(mapply(function(featureMeasurements, scale, location)
  sum(dnorm(featureMeasurements, scale, location, log = TRUE)),
  measurements, allDistribution[[1]], allDistribution[[2]])) -
  rowSums(sapply(levels(classes), function(class)
  {
    classMeasurements <- measurements[which(classes == class), ]
    classDistribution <- getLocationsAndScales(classMeasurements, ...)
    
    unlist(mapply(function(featureMeasurements, scale, location)
    sum(dnorm(featureMeasurements, scale, location, log = TRUE)),
    classMeasurements,
    switch(alternative[["location"]], same = allDistribution[[1]], different = classDistribution[[1]]),
    switch(alternative[["scale"]], same = allDistribution[[2]], different = classDistribution[[2]])))    
  }))
  
  if(!is.null(S4Vectors::mcols(measurements)))
    S4Vectors::mcols(measurements)[order(logLikelihoodRatios), ]
  else
    colnames(measurements)[order(logLikelihoodRatios)]
})

# One or more omics data sets, possibly with clinical data.
setMethod("likelihoodRatioRanking", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), classes, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]

  if(ncol(dataTable) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
  likelihoodRatioRanking(dataTable, classes, ...)
})