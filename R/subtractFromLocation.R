#' Subtract Numeric Feature Measurements from a Location
#' 
#' For each numeric feature, calculates the location, and subtracts all
#' measurements from that location.
#' 
#' Only the samples specified by \code{measurementsTrain} are used in the calculation of
#' the location.
#' 
#' @aliases subtractFromLocation subtractFromLocation,matrix,matrix-method
#' subtractFromLocation,DataFrame,DataFrame-method
#' subtractFromLocation,MultiAssayExperiment,MultiAssayExperiment-method
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
#' If of type \code{\link{DataFrame}} or \code{\link{MultiAssayExperiment}}, the data set is subset
#' to only those features of type \code{numeric}.
#' @param measurementsTest A data set of the same type as \code{measurementsTrain} with no samples in common with it.
#' The subtraction will also be performed to it.
#' @param location Character. Either "mean" or "median".
#' @param absolute Logical. Default: \code{TRUE}. If \code{TRUE}, then absolute
#' values of the differences are returned. Otherwise, they are signed.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"sampleInfo"} is also a valid value
#' and specifies that numeric variables from the sample information data table will be
#' used.
#' @param verbose Default: 3. A progress message is shown if this value is 3.
#' @return The same class of variable as the input variable \code{measurements}
#' is, with the numeric features subtracted from the calculated location.
#' @author Dario Strbenac
#' @examples
#'   aMatrix <- matrix(1:100, ncol = 10)
#'   subtractFromLocation(aMatrix[1:5,], aMatrix[6:10, ], "median")
#' @export
#' @usage NULL
setGeneric("subtractFromLocation", function(measurementsTrain, measurementsTest, ...) standardGeneric("subtractFromLocation"))

#' @export
setMethod("subtractFromLocation", c("matrix", "matrix"), function(measurementsTrain, measurementsTest, location = c("mean", "median"),
           absolute = TRUE, verbose = 3)
{
  location <- match.arg(location)
  if(location == "mean")
    featureTrainingLocations <- colMeans(measurementsTrain, na.rm = TRUE)
  else # median.
    featureTrainingLocations <- apply(measurementsTrain, 2, median, na.rm = TRUE)
  
  transformedTrain <- t(apply(measurementsTrain, 1, '-', featureTrainingLocations))
  transformedTest <- t(apply(measurementsTest, 1, '-', featureTrainingLocations))
  if(absolute == TRUE)
  {
    transformedTrain <- abs(transformedTrain)
    transformedTest <- abs(transformedTest)
  }
  
  if(verbose == 3)
    message("Subtraction from ", location,
            {if(absolute == TRUE) " and absolute transformation"}, " completed.")
            
  list(transformedTrain, transformedTest)
})

# Sample information data or one of the other inputs, transformed.
#' @export
setMethod("subtractFromLocation", c("DataFrame", "DataFrame"), function(measurementsTrain, measurementsTest, location = c("mean", "median"),
           absolute = TRUE, verbose = 3)
{
  isNumeric <- sapply(measurementsTrain, is.numeric)
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  
  if(sum(isNumeric) != ncol(measurementsTrain) && verbose == 3)
    message("Some columns are not numeric. Only subtracting values from location for\n", 
            "the columns containing numeric data.")
  
  location <- match.arg(location)
  measurementsTrain <- measurementsTrain[, isNumeric]
  measurementsTest <- measurementsTest[, isNumeric]
  if(location == "mean")
    locations <- apply(measurementsTrain, 2, mean, na.rm = TRUE)
  else # median.
    locations <- apply(measurementsTrain, 2, median, na.rm = TRUE)
  
  transformedTrain <- DataFrame(t(apply(measurementsTrain, 1, '-', locations)))
  transformedTest <- DataFrame(t(apply(measurementsTest, 1, '-', locations)))
  if(absolute == TRUE)
  {
    transformedTrain <- DataFrame(lapply(transformedTrain, abs))
    transformedTest <- DataFrame(lapply(transformedTest, abs))
  }
  
  if(verbose == 3)
    message("Subtraction from ", location,
            {if(absolute == TRUE) " and absolute transformation"}, " completed.")
  
  list(transformedTrain, transformedTest)
})

#' @export
setMethod("subtractFromLocation", "MultiAssayExperiment",
          function(measurementsTrain, measurementsTest, targets = names(measurementsTrain),
           location = c("mean", "median"), absolute = TRUE, verbose = 3)
{
  location <- match.arg(location)
  if(!all(targets %in% c(names(measurementsTrain), "sampleInfo")))
    stop("Some table names in 'targets' are not assay names in 'measurementsTrain' or \"sampleInfo\".")  
  
  if("sampleInfo" %in% targets)
  {
    targets <- targets[-which(match, "sampleInfo")]
    isNumeric <- sapply(MultiAssayExperiment::colData(measurementsTrain), is.numeric)
    sampleInfoTrain <- MultiAssayExperiment::colData(measurementsTrain)[, isNumeric]
    if(location == "mean")
      locations <- apply(sampleInfoTrain, 2, mean, na.rm = TRUE)
    else # median.
      locations <- apply(sampleInfoTrain, 2, median, na.rm = TRUE)
    transformedSampleInfoTrain <- DataFrame(t(apply(MultiAssayExperiment::colData(measurementsTrain), 1, '-', locations)))
    transformedSampleInfoTest <- DataFrame(t(apply(MultiAssayExperiment::colData(measurementsTest), 1, '-', locations)))
    if(absolute == TRUE)
    {
      transformedSampleInfoTrain <- DataFrame(lapply(transformedSampleInfoTrain, abs))
      transformedSampleInfoTest <- DataFrame(lapply(transformedSampleInfoTest, abs))
    }
    MultiAssayExperiment::colData(measurementsTrain)[, isNumeric] <- transformedSampleInfoTrain
    MultiAssayExperiment::colData(measurementsTest)[, isNumeric] <- transformedSampleInfoTest
  }
  
  # Now, transform the experiments. The features are in the rows, unlike for colData.
  measurementsTrain <- measurementsTrain[, , targets]
  if(location == "mean")
    featureTrainingLocations <- lapply(MultiAssayExperiment::experiments(measurementsTrain), rowMeans, na.rm = TRUE)
  else # median.
    featureTrainingLocations <- lapply(MultiAssayExperiment::experiments(measurementsTrain), function(dataTable)
                                                apply(dataTable, 1, median, na.rm = TRUE))
  transformedTablesTrain <- mapply(function(dataTable, locations) apply(dataTable, 2, '-', locations),
              MultiAssayExperiment::experiments(measurementsTrain[, , targets]), featureTrainingLocations, SIMPLIFY = FALSE)
  transformedTablesTest <- mapply(function(dataTable, locations) apply(dataTable, 2, '-', locations),
              MultiAssayExperiment::experiments(measurementsTest[, , targets]), featureTrainingLocations, SIMPLIFY = FALSE)
  if(absolute == TRUE)
  {
    transformedTablesTrain <- lapply(transformedTablesTrain, abs)
    transformedTablesTest <- lapply(transformedTablesTest, abs)
  }
  
  MultiAssayExperiment::experiments(measurementsTrain)[targets] <- MultiAssayExperiment::ExperimentList(transformedTablesTrain)
  MultiAssayExperiment::experiments(measurementsTest)[targets] <- MultiAssayExperiment::ExperimentList(transformedTablesTest)

  if(verbose == 3)
    message("Subtraction from ", location,
            {if(absolute == TRUE) " and absolute transformation"}, " completed.")
  
  list(measurementsTrain, measurementsTest)
})