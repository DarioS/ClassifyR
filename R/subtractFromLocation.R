#' Subtract Numeric Feature Measurements from a Location
#' 
#' For each numeric feature, calculates the location, and subtracts all
#' measurements from that location.
#' 
#' Only the samples specified by \code{training} are used in the calculation of
#' the location.  To use all samples for calculation of the location, simply
#' provide indices of all the samples.
#' 
#' @aliases subtractFromLocation subtractFromLocation,matrix-method
#' subtractFromLocation,DataFrame-method
#' subtractFromLocation,MultiAssayExperiment-method
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
#' If of type \code{\link{DataFrame}} or \code{\link{MultiAssayExperiment}}, the data set is subset
#' to only those features of type \code{numeric}.
#' @param location Character. Either "mean" or "median".
#' @param absolute Logical. Default: \code{TRUE}. If \code{TRUE}, then absolute
#' values of the differences are returned. Otherwise, they are signed.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"sampleInfo"} is also a valid value
#' and specifies that numeric variables from the sample information data table will be
#' used.
#' @param transformName Default: Location Subtraction. Useful for automated
#' plot annotation by plotting functions within this package.
#' @param verbose Default: 3. A progress message is shown if this value is 3.
#' @return The same class of variable as the input variable \code{measurements}
#' is, with the numeric features subtracted from the calculated location.
#' @author Dario Strbenac
#' @examples
#' 
#'   aMatrix <- matrix(1:100, ncol = 10)
#'   subtractFromLocation(aMatrix, training = 1:5, "median")
#' 
#' @importFrom MultiAssayExperiment ExperimentList colData experiments MultiAssayExperiment
#' @export
setGeneric("subtractFromLocation", function(measurementsTrain, ...) standardGeneric("subtractFromLocation"))

setMethod("subtractFromLocation", "matrix", function(measurementsTrain, location = c("mean", "median"),
           absolute = TRUE, transformName = "Location Subtraction", verbose = 3)
{
  location <- match.arg(location)
  if(location == "mean")
    featureTrainingLocations <- colMeans(measurementsTrain, na.rm = TRUE)
  else # median.
    featureTrainingLocations <- apply(measurementsTrain, 2, median, na.rm = TRUE)
  transformed <- apply(measurementsTrain, 1, '-', featureTrainingLocations)
  if(absolute == TRUE)
    transformed <- abs(transformed)
  
  if(verbose == 3)
    message("Subtraction from ", location,
            {if(absolute == TRUE) " and absolute transformation"}, " completed.")
            
  transformed # Return an endomorphic variable; a matrix.
})

# Sample information data or one of the other inputs, transformed.
setMethod("subtractFromLocation", "DataFrame", function(measurementsTrain, location = c("mean", "median"),
           absolute = TRUE, transformName = "Location Subtraction", verbose = 3)
{
  isNumeric <- sapply(measurementsTrain, is.numeric)
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  
  if(sum(isNumeric) != ncol(measurementsTrain) && verbose == 3)
    message("Some columns are not numeric. Only subtracting values from location for\n", 
            "the columns containing numeric data.")
  
  location <- match.arg(location)
  measurementsTrain <- measurementsTrain[, isNumeric]
  if(location == "mean")
    locations <- apply(measurementsTrain, 2, mean, na.rm = TRUE)
  else # median.
    locations <- apply(measurementsTrain, 2, median, na.rm = TRUE)
  transformed <- measurementsTrain
  transformed[, isNumeric] <- DataFrame(t(apply(measurementsTrain[, isNumeric], 1, '-', locations)))
  if(absolute == TRUE)
    transformed[, isNumeric] <- DataFrame(lapply(transformed[, isNumeric], abs))
  
  if(verbose == 3)
    message("Subtraction from ", location,
            {if(absolute == TRUE) " and absolute transformation"}, " completed.")
  
  transformed # Return an endomorphic variable; a DataFrame.
})

setMethod("subtractFromLocation", "MultiAssayExperiment", function(measurementsTrain, targets = names(measurementsTrain),
           location = c("mean", "median"), absolute = TRUE, transformName = "Location Subtraction", verbose = 3)
{
  location <- match.arg(location)
  if(!all(targets %in% c(names(measurementsTrain), "sampleInfo")))
    stop("Some table names in 'targets' are not assay names in 'measurementsTrain' or \"sampleInfo\".")  
  
  transformed <- measurementsTrain
  
  if("sampleInfo" %in% targets)
  {
    targets <- targets[-which(match, "sampleInfo")]
    isNumeric <- sapply(MultiAssayExperiment::colData(measurementsTrain), is.numeric)
    sampleInfoTrain <- MultiAssayExperiment::colData(measurementsTrain)[, isNumeric]
    if(location == "mean")
      locations <- apply(sampleInfoTrain, 2, mean, na.rm = TRUE)
    else # median.
      locations <- apply(sampleInfoTrain, 2, median, na.rm = TRUE)
    transformedSampleInfo <- DataFrame(t(apply(MultiAssayExperiment::colData(measurementsTrain), 1, '-', locations)))
    if(absolute == TRUE)
      transformedSampleInfo <- DataFrame(lapply(transformedSampleInfo, abs))
    MultiAssayExperiment::colData(transformed)[, isNumeric] <- transformedSampleInfo
  }
  
  measurementsTrain <- measurementsTrain[, , targets]
  if(location == "mean")
    featureTrainingLocations <- lapply(MultiAssayExperiment::experiments(measurementsTrain), rowMeans, na.rm = TRUE)
  else # median.
    featureTrainingLocations <- lapply(MultiAssayExperiment::experiments(measurementsTrain), function(dataTable)
                                                apply(dataTable, 1, median, na.rm = TRUE))
  transformedTables <- mapply(function(dataTable, locations) apply(dataTable, 2, '-', locations),
              MultiAssayExperiment::experiments(measurementsTrain[, , targets]), featureTrainingLocations, SIMPLIFY = FALSE)
  if(absolute == TRUE)
    transformedTables <- lapply(transformedTables, abs)
  
  MultiAssayExperiment::experiments(transformed)[targets] <- MultiAssayExperiment::ExperimentList(transformedTables)

  if(verbose == 3)
    message("Subtraction from ", location,
            {if(absolute == TRUE) " and absolute transformation"}, " completed.")
  
  transformed # Return an endomorphic variable; a MultiAssayExperiment.
})