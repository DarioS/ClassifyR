################################################################################
#
# Train interface
#
################################################################################


#' An Interface for randomForestSRC Package's rfsrc random forest survival Function
#' 
#' Random Forest Survival Model
#' 
#' @aliases rfsrcInterface rfsrcTrainInterface rfsrcPredictInterface
#' rfsrcInterface,matrix-method
#' rfsrcInterface,DataFrame-method
#' rfsrcInterface,MultiAssayExperiment-method
#' rfsrcPredictInterface
#' rfsrcPredictInterface,rfsrc,matrix-method
#' rfsrcPredictInterface,rfsrc,DataFrame-method
#' rfsrcPredictInterface,rfsrc,MultiAssayExperiment-method
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
#' @param survivalTrain A tabular data type of survival information of the
#' same number of rows as the number of samples in \code{measurementsTrain} and 2 to 3 columns if it is a
#' \code{\link{matrix}} or a \code{\link{DataFrame}}, or a character vector of length 2 to 3 containing the
#' column names in \code{measurementsTrain} if it is a \code{\link{DataFrame}} or the
#' column name in \code{colData(measurementsTrain)} if \code{measurementsTrain} is a
#' \code{\link{MultiAssayExperiment}}. If a vector of column names, those columns will be
#' removed before training.
#' @param model A trained rfsrc classifier, as created by
#' \code{rfsrcTrainInterface}, which has the same form as the output of
#' \code{\link[randomForestSRC]{rfsrc}}.
#' @param measurementsTest An object of the same class as \code{measurementsTrain} with no
#' samples in common with \code{measurementsTrain} and the same number of features
#' as it.
#' @param targets If \code{measurementsTrain} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"sampleInfo"} is also a valid value
#' and specifies that integer variables from the clinical data table will be
#' used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method (e.g. \code{verbose}) or options which are accepted
#' by the \code{\link[randomForestSRC]{rfsrc}} or \code{\link[randomForestSRC]{predict.rfsrc}} functions.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return For \code{rfsrcTrainInterface}, the trained Random Forest Survival model.
#' For \code{rfsrcPredictInterface}, a risk score prediction for each sample.
#' 
#' @importFrom survival Surv concordance
#' @rdname rfsrcInterface
#' @usage NULL
#' @export
setGeneric("rfsrcTrainInterface", function(measurementsTrain, ...) standardGeneric("rfsrcTrainInterface"))

#' @rdname rfsrcInterface
#' @export
setMethod("rfsrcTrainInterface", "matrix", function(measurementsTrain, survivalTrain, ...)
{
  rfsrcTrainInterface(S4Vectors::DataFrame(measurementsTrain, check.names = FALSE), survivalTrain, ...)
})

# Clinical data or one of the other inputs, transformed.
#' @rdname rfsrcInterface
#' @export
setMethod("rfsrcTrainInterface", "DataFrame", function(measurementsTrain, survivalTrain, ..., verbose = 3)
{
  if(!requireNamespace("survival", quietly = TRUE))
    stop("The package 'survival' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting coxph classifier to training data and making predictions on test
            data.")
  
  if(!requireNamespace("randomForestSRC", quietly = TRUE))
    stop("The package 'randomForestSRC' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting rfsrc classifier to training data and making predictions on test
            data.")

  splitDataset <- .splitDataAndOutcome(measurementsTrain, survivalTrain)  
  survivalTrain <- splitDataset[["outcome"]]
  measurementsTrain <- splitDataset[["measurements"]]
  bindedMeasurements <- cbind(measurementsTrain, event = survivalTrain[,1], time = survivalTrain[,2])
  randomForestSRC::rfsrc(Surv(event = event, time = time) ~ ., as.data.frame(bindedMeasurements), ...)
})

#' @rdname rfsrcInterface
#' @export
setMethod("rfsrcTrainInterface", "MultiAssayExperiment", function(measurementsTrain, targets = names(measurementsTrain), survivalTrain, ...)
{
  tablesAndSurvival <- ClassifyR:::.MAEtoWideTable(measurementsTrain, targets, survivalTrain, restrict = NULL)
  measurementsTrain <- tablesAndSurvival[["dataTable"]]
  survivalTrain <- tablesAndSurvival[["outcome"]]
  
  rfsrcTrainInterface(measurementsTrain, survivalTrain, ...)
})



################################################################################
#
# Predict Interface
#
################################################################################


#' @rdname rfsrcInterface
#' @usage NULL
#' @export
setGeneric("rfsrcPredictInterface", function(model, measurementsTest, ...)
  standardGeneric("rfsrcPredictInterface"))

#' @rdname rfsrcInterface
#' @export
setMethod("rfsrcPredictInterface", c("rfsrc", "matrix"), # Matrix of numeric measurements.
          function(model, measurementsTest, ...)
{
  rfsrcPredictInterface(model, S4Vectors::DataFrame(measurementsTest, check.names = FALSE), ...)
})

#' @rdname rfsrcInterface
#' @export
setMethod("rfsrcPredictInterface", c("rfsrc", "DataFrame"),
function(model, measurementsTest, ..., verbose = 3)
{
  predictedOutcome = predict(model, as.data.frame(measurementsTest), ...)$predicted
  names(predictedOutcome) = rownames(measurementsTest)
  predictedOutcome
})

# One or more omics data sets, possibly with clinical data.
#' @rdname rfsrcInterface
#' @export
setMethod("rfsrcPredictInterface", c("rfsrc", "MultiAssayExperiment"),
          function(model, measurementsTest, targets = names(measurementsTest), ...)
{
  testingTable <- ClassifyR:::.MAEtoWideTable(measurementsTest, targets)
  rfsrcPredictInterface(model, testingTable, ...)
})