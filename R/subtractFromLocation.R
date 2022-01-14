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
#' @param measurements A \code{\link{matrix}}, \code{\link{DataFrame}} or a
#' \code{\link{MultiAssayExperiment}} object containing the data.  For a
#' matrix, the rows are for features and the columns are for samples.
#' @param training A vector specifying which samples are in the training set.
#' @param location Character. Either "mean" or "median".
#' @param absolute Logical. Default: \code{TRUE}. If \code{TRUE}, then absolute
#' values of the differences are returned. Otherwise, they are signed.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"clinical"} is also a valid value
#' and specifies that numeric variables from the clinical data table will be
#' used.
#' @param transformName Default: Location Subtraction. Useful for automated
#' plot annotation by plotting functions within this package..
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
setGeneric("subtractFromLocation", function(measurements, ...)
           standardGeneric("subtractFromLocation"))

setMethod("subtractFromLocation", "matrix", 
          function(measurements, training, location = c("mean", "median"),
                   absolute = TRUE, transformName = "Location Subtraction", verbose = 3)
{
  location <- match.arg(location)
  measurementsTrain <- measurements[, training]
  if(location == "mean")
    featureTrainingLocations <- rowMeans(measurementsTrain, na.rm = TRUE)
  else # median.
    featureTrainingLocations <- apply(measurementsTrain, 1, median, na.rm = TRUE)
  transformed <- apply(measurements, 2, '-', featureTrainingLocations)
  if(absolute == TRUE)
    transformed <- abs(transformed)
  
  if(verbose == 3)
    message("Subtraction from ", location,
            {if(absolute == TRUE) " and absolute transformation"}, " completed.")
            
  transformed # Return an endomorphic variable; a matrix.
})

setMethod("subtractFromLocation", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, training, location = c("mean", "median"),
                   absolute = TRUE, transformName = "Location Subtraction", verbose = 3)
          {
            isNumeric <- sapply(measurements, is.numeric)
            if(sum(isNumeric) == 0)
              stop("No features are numeric but at least one must be.")
            
            if(sum(isNumeric) != ncol(measurements) && verbose == 3)
              message("Some columns are not numeric. Only subtracting values from location for\n", 
                      "the columns containing numeric data.")
            
            location <- match.arg(location)
            measurementsTrain <- measurements[training, isNumeric]
            if(location == "mean")
              locations <- apply(measurementsTrain, 2, mean, na.rm = TRUE)
            else # median.
              locations <- apply(measurementsTrain, 2, median, na.rm = TRUE)
            transformed <- measurements
            transformed[, isNumeric] <- DataFrame(t(apply(measurements[, isNumeric], 1, '-', locations)))
            if(absolute == TRUE)
              transformed[, isNumeric] <- DataFrame(lapply(transformed[, isNumeric], abs))
            
            if(verbose == 3)
              message("Subtraction from ", location,
                      {if(absolute == TRUE) " and absolute transformation"}, " completed.")
            
            transformed # Return an endomorphic variable; a DataFrame.
          })

setMethod("subtractFromLocation", "MultiAssayExperiment", 
          function(measurements, training, targets = names(measurements),
                   location = c("mean", "median"), absolute = TRUE, transformName = "Location Subtraction", verbose = 3)
{
  location <- match.arg(location)
  if(!all(targets %in% c(names(measurements), "clinical")))
    stop("Some table names in 'targets' are not assay names in 'measurements' or \"clinical\".")  
  
  transformed <- measurements
  
  if("clinical" %in% targets)
  {
    targets <- targets[-which(match, "clinical")]
    isNumeric <- sapply(MultiAssayExperiment::colData(measurements), is.numeric)
    clinicalTrain <- MultiAssayExperiment::colData(measurements)[training, isNumeric]
    if(location == "mean")
      locations <- apply(clinicalTrain, 2, mean, na.rm = TRUE)
    else # median.
      locations <- apply(clinicalTrain, 2, median, na.rm = TRUE)
    transformedClinical <- DataFrame(t(apply(MultiAssayExperiment::colData(measurements), 1, '-', locations)))
    if(absolute == TRUE)
      transformedClinical <- DataFrame(lapply(transformedClinical, abs))
    MultiAssayExperiment::colData(transformed)[, isNumeric] <- transformedClinical
  }
  
  measurementsTrain <- measurements[, training, targets]
  if(location == "mean")
    featureTrainingLocations <- lapply(MultiAssayExperiment::experiments(measurementsTrain), rowMeans, na.rm = TRUE)
  else # median.
    featureTrainingLocations <- lapply(MultiAssayExperiment::experiments(measurementsTrain), function(dataTable)
                                                apply(dataTable, 1, median, na.rm = TRUE))
  transformedTables <- mapply(function(dataTable, locations) apply(dataTable, 2, '-', locations),
              MultiAssayExperiment::experiments(measurements[, , targets]), featureTrainingLocations, SIMPLIFY = FALSE)
  if(absolute == TRUE)
    transformedTables <- lapply(transformedTables, abs)
  
  MultiAssayExperiment::experiments(transformed)[targets] <- MultiAssayExperiment::ExperimentList(transformedTables)

  if(verbose == 3)
    message("Subtraction from ", location,
            {if(absolute == TRUE) " and absolute transformation"}, " completed.")
  
  transformed # Return an endomorphic variable; a MultiAssayExperiment.
})