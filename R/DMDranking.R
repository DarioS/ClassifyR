
# Mean or Median, Three Kinds of Deviations.

#' Ranking by Differential Distributions with Differences in Means or Medians
#' and a Deviation Measure
#' 
#' Ranks features by largest Differences in Means/Medians and Deviations.
#' 
#' DMD is defined as \eqn{\sum_{i = 1}\sum_{j = i + 1}|location_i - location_j|
#' + |scale_i - scale_j|}{sum(|location i - location j|) + sum(|scale i - scale
#' j|), i not j, i < j}. The subscripts denote the class for which the
#' parameter is calculated for.
#' 
#' @aliases DMDranking DMDranking,matrix-method DMDranking,DataFrame-method
#' DMDranking,MultiAssayExperiment-method
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix}, the rows are features, and the columns are samples.
#' @param classes A vector of class labels of class \code{\link{factor}} of the
#' same length as the number of samples in \code{measurements} if it is a
#' \code{\link{matrix}} (i.e. number of columns) or a \code{\link{DataFrame}}
#' (i.e. number of rows) or a character vector of length 1 containing the
#' column name in \code{measurements} if it is a \code{\link{DataFrame}} or the
#' column name in \code{colData(measurements)} if \code{measurements} is a
#' \code{\link{MultiAssayExperiment}}. If a column name, that column will be
#' removed before training.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"clinical"} is also a valid value
#' and specifies that numeric variables from the clinical data table will be
#' used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method or parameters for
#' \code{\link{getLocationsAndScales}}, such as \code{location}, \code{scale}.
#' @param differences Default: \code{"both"}. Either \code{"both"},
#' \code{"location"}, or \code{"scale"}.  The type of differences to consider.
#' If both are considered then the absolute difference in location and the
#' absolute difference in scale are summed.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return A vector of feature indices, from the most promising features in the
#' first position to the least promising feature in the last position.
#' @author Dario Strbenac
#' @examples
#' 
#'   # First 20 features have bimodal distribution for Poor class.
#'   # Other 80 features have normal distribution for both classes.
#'   genesMatrix <- sapply(1:25, function(sample)
#'                               {
#'                                 randomMeans <- sample(c(8, 12), 20, replace = TRUE)
#'                                 c(rnorm(20, randomMeans, 1), rnorm(80, 10, 1))
#'                               }
#'                        )
#'   genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample) rnorm(100, 10, 1)))
#'   classes <- factor(rep(c("Poor", "Good"), each = 25))
#'   
#'   ranked <- DMDranking(genesMatrix, classes)
#'   head(ranked)
#' 
#' @export
setGeneric("DMDranking", function(measurements, ...)
           standardGeneric("DMDranking"))

setMethod("DMDranking", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  DMDranking(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("DMDranking", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, differences = c("both", "location", "scale"),
                   ..., verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  
  if(verbose == 3)
    message("Selecting features by DMD.")
  differences <- match.arg(differences)
  
  allClassesLocationsScales <- lapply(levels(classes), function(class)
  {
    aClassMeasurements <- measurements[which(classes == class), ]
    getLocationsAndScales(aClassMeasurements, ...)
  })
  allClassesLocations <- sapply(allClassesLocationsScales, "[[", 1)
  allClassesScales <- sapply(allClassesLocationsScales, "[[", 2)
  locationsDifferences <- apply(allClassesLocations, 1, function(locations) sum(abs(c(dist(locations)))))
  scalesDifferences <- apply(allClassesScales, 1, function(scales) sum(abs(c(dist(scales)))))
  
  divergence <- 0
  if(differences %in% c("both", "location"))
    divergence <- divergence + locationsDifferences
  if(differences %in% c("both", "scale"))
    divergence <- divergence + scalesDifferences

  order(divergence, decreasing = TRUE)
})

# One or more omics data sets, possibly with clinical data.
setMethod("DMDranking", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), classes, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]            
  DMDranking(dataTable, classes, ...)
})