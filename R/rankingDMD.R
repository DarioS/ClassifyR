
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
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
#' If of type \code{\link{DataFrame}} or \code{\link{MultiAssayExperiment}}, the data set is subset
#' to only those features of type \code{numeric}.
#' @param classesTrain A vector of class labels of class \code{\link{factor}} of the
#' same length as the number of samples in \code{measurementsTrain} if it is a
#' \code{\link{matrix}} or a \code{\link{DataFrame}} or a character vector of length 1
#' containing the column name in \code{measurementsTrain} if it is a \code{\link{DataFrame}} or the
#' column name in \code{colData(measurementsTrain)} if \code{measurementsTrain} is a
#' \code{\link{MultiAssayExperiment}}. If a column name, that column will be
#' removed before training.
#' @param targets If \code{measurementsTrain} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"sampleInfo"} is also a valid value
#' and specifies that numeric variables from the sample information data table will be
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
#'   genesMatrix <- sapply(1:20, function(feature)
#'                               {
#'                                 randomMeans <- sample(c(8, 12), 25, replace = TRUE)
#'                                 c(rnorm(25, randomMeans, 1), rnorm(25, 10, 1))
#'                               }
#'                        )
#'   genesMatrix <- cbind(genesMatrix, sapply(1:80, function(feature) rnorm(50, 10, 1)))
#'   classes <- factor(rep(c("Poor", "Good"), each = 25))
#'   
#'   ranked <- DMDranking(genesMatrix, classes)
#'   head(ranked)
#' 
#' @rdname DMDranking
#' @usage NULL
#' @export
setGeneric("DMDranking", function(measurementsTrain, ...)
           standardGeneric("DMDranking"))

#' @rdname DMDranking
#' @export
setMethod("DMDranking", "matrix", # Matrix of numeric measurements.
          function(measurementsTrain, classesTrain, ...)
{
  DMDranking(DataFrame(measurementsTrain, check.names = FALSE), classesTrain, ...)
})

#' @rdname DMDranking
#' @export
setMethod("DMDranking", "DataFrame", # sampleInfo data or one of the other inputs, transformed.
          function(measurementsTrain, classesTrain, differences = c("both", "location", "scale"),
                   ..., verbose = 3)
{
  splitDataset <- .splitDataAndOutcomes(measurementsTrain, classesTrain)
  measurementsTrain <- splitDataset[["measurements"]]

  if(verbose == 3)
    message("Selecting features by DMD.")
  differences <- match.arg(differences)
  
  allClassesLocationsScales <- lapply(levels(classesTrain), function(class)
  {
    aClassMeasurements <- measurementsTrain[which(classesTrain == class), ]
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

#' @rdname DMDranking
#' @export
# One or more omics data sets, possibly with sample information data.
setMethod("DMDranking", "MultiAssayExperiment",
          function(measurementsTrain, targets = names(measurementsTrain), classesTrain, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurementsTrain, targets, classesTrain)
  measurementsTrain <- tablesAndClasses[["dataTable"]]
  classesTrain <- tablesAndClasses[["outcomes"]]            
  DMDranking(measurementsTrain, classesTrain, ...)
})