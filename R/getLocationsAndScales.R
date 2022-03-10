#' Calculate Location and Scale
#' 
#' Calculates the location and scale for each feature.
#' 
#' \code{"SD"} is used to represent standard deviation and \code{"MAD"} is used
#' to represent median absolute deviation.
#' 
#' @aliases getLocationsAndScales getLocationsAndScales,matrix-method
#' getLocationsAndScales,DataFrame-method
#' getLocationsAndScales,MultiAssayExperiment-method
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
#' \code{matrix}, the rows are samples, and the columns are features.
#' If of type \code{\link{DataFrame}} or \code{\link{MultiAssayExperiment}}, the data set is subset
#' to only those features of type \code{numeric}.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"sampleInfo"} is also a valid value
#' and specifies that numeric variables from the sample information data table will be
#' used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method.
#' @param location The type of location to be calculated.
#' @param scale The type of scale to be calculated.
#' @return A \code{\link{list}} of length 2. The first element contains the
#' location for every feature. The second element contains the scale for every
#' feature.
#' @author Dario Strbenac
#' @references Qn:
#' \url{http://www.tandfonline.com/doi/pdf/10.1080/01621459.1993.10476408}
#' @examples
#' 
#'   genesMatrix <- matrix(rnorm(1000, 8, 4), nrow = 10)
#'   distributionInfo <- getLocationsAndScales(genesMatrix, "median", "MAD")
#'   
#'   mean(distributionInfo[["median"]]) # Typical median.
#'   mean(distributionInfo[["MAD"]]) # Typical MAD.
#' 
#' @usage NULL
#' @export
setGeneric("getLocationsAndScales", function(measurements, ...)
           standardGeneric("getLocationsAndScales"))

#' @export
setMethod("getLocationsAndScales", "matrix", # Matrix of numeric measurements.
          function(measurements, ...)
{
  getLocationsAndScales(DataFrame(measurements, check.names = FALSE), ...)
})

#' @export
setMethod("getLocationsAndScales", "DataFrame", # Sample information data or one of the other inputs, transformed.
          function(measurements, location = c("mean", "median"), scale = c("SD", "MAD", "Qn"))
{
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  location <- match.arg(location) 
  scale <- match.arg(scale)
   
  if(scale == "Qn" && !requireNamespace("robustbase", quietly = TRUE))
    stop("The package 'robustbase' could not be found. Please install it.")
 
  setNames(list(switch(location,
                       mean = apply(measurements, 2, mean),
                       median = apply(measurements, 2, median)),
                switch(scale,
                       SD = apply(measurements, 2, sd),
                       MAD = apply(measurements, 2, mad),
                       Qn = apply(measurements, 2, robustbase::Qn))),
                c(location, scale))
})

# One or more omics data sets, possibly with sample information data.
#' @export
setMethod("getLocationsAndScales", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
{
  if(!all(targets %in% c(names(measurements), "sampleInfo")))
    stop("Some table names in 'targets' are not assay names in 'measurements' or \"sampleInfo\".")  
            
  combinedData <- .MAEtoWideTable(measurements, targets, NULL)
  if(class(combinedData) == "list")
    combinedData <- combinedData[["dataTable"]]
  getLocationsAndScales(combinedData, ...)            
})