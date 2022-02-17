################################################################################
#
# Train Interface
#
################################################################################




#' An Interface for glmnet Package's coxnet Function
#' 
#' An elastic net GLM classifier uses a penalty which is a combination of a
#' lasso penalty and a ridge penalty, scaled by a lambda value, to fit a sparse
#' linear model to the data.
#' 
#' The value of the \code{family} parameter is fixed to \code{"cox"} so
#' that classification with survival is possible.
#' During classifier training, if more than one lambda value
#' is considered by specifying a vector of them as input or leaving the default
#' value of NULL, then the chosen value is determined based on classifier
#' resubstitution error rate.
#' 
#' @aliases coxnetInterface coxnetTrainInterface
#' coxnetPredictInterface coxnetTrainInterface,matrix-method
#' coxnetTrainInterface,DataFrame-method
#' coxnetTrainInterface,MultiAssayExperiment-method
#' coxnetPredictInterface,multnet,matrix-method
#' coxnetPredictInterface,multnet,DataFrame-method
#' coxnetPredictInterface,multnet,MultiAssayExperiment-method
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
#' @param survivalTrain A tabular data type of survival information of the
#' same number of rows as the number of samples in \code{measurementsTrain} and 2 to 3 columns if it is a
#' \code{\link{matrix}} or a \code{\link{DataFrame}}, or a character vector of length 2 to 3 containing the
#' column names in \code{measurementsTrain} if it is a \code{\link{DataFrame}} or the
#' column name in \code{colData(measurementsTrain)} if \code{measurementsTrain} is a
#' \code{\link{MultiAssayExperiment}}. If a vector of column names, those columns will be
#' removed before training.
#' @param lambda The lambda value passed directly to
#' \code{\link[glmnet]{glmnet}} if the training function is used or passed as
#' \code{s} to \code{\link[glmnet]{predict.glmnet}} if the prediction function
#' is used.
#' @param measurementsTest An object of the same class as \code{measurementsTrain} with no
#' samples in common with \code{measurements} and the same number of features
#' as it.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"clinical"} is also a valid value
#' and specifies that integer variables from the clinical data table will be
#' used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method (e.g. \code{verbose}) or, for the training function,
#' options that are used by the \code{glmnet} function. For the testing
#' function, this variable simply contains any parameters passed from the
#' classification framework to it which aren't used by glmnet's \code{predict}
#' fuction.
#' @param model A trained coxnet, as created by the \code{glmnet}
#' function.
#' @param returnType Default: \code{"both"}. Either \code{"class"},
#' \code{"score"} or \code{"both"}.  Sets the return value from the prediction
#' to either a vector of class labels, matrix of scores for each class, or both
#' labels and scores in a \code{data.frame}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return For \code{coxnetTrainInterface}, an object of type
#' \code{glmnet}. For \code{coxnetPredictInterface}, a matrix of containing the link and risk functions.
#' \code{returnType}.
#' @examples
#' 
#'   
#' @export
setGeneric("coxnetTrainInterface", function(measurementsTrain, ...)
  standardGeneric("coxnetTrainInterface"))

setMethod("coxnetTrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurementsTrain, survivalTrain, ...)
          {
            coxnetTrainInterface(DataFrame(measurementsTrain, check.names = FALSE), survivalTrain, ...)
          })

# Clinical data or one of the other inputs, transformed.
setMethod("coxnetTrainInterface", "DataFrame", function(measurementsTrain, survivalTrain, lambda = NULL, ..., verbose = 3)
{
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting coxnet model to data.")
  
  splitDataset <- .splitDataAndOutcomes(measurementsTrain, survivalTrain)
  measurementsTrain <- data.frame(splitDataset[["measurements"]], check.names = FALSE)
  measurementsMatrix <- glmnet::makeX(as(measurementsTrain, "data.frame"))
  
  # The response variable is a Surv class of object.
  fit <- glmnet::cv.glmnet(measurementsMatrix, splitDataset[["outcomes"]], family = "cox", type = "C", ...)
  fitted <- fit$glmnet.fit
  
  offset <- -mean(predict(fitted, measurementsMatrix, s = fit$lambda.min, type = "link"))
  attr(fitted, "tune") <- list(lambda = fit$lambda.min, offset = offset)
  
  fitted
})

# One or more omics datasets, possibly with clinical data.
setMethod("coxnetTrainInterface", "MultiAssayExperiment",
          function(measurementsTrain, targets = names(measurementsTrain), survivalTrain, ...)
          {
            tablesAndClasses <- .MAEtoWideTable(measurementsTrain, targets, survivalTrain)
            measurementsTrain <- tablesAndClasses[["dataTable"]]
            survivalTrain <- tablesAndClasses[["outcomes"]]
            
            if(ncol(measurementsTrain) == 0)
              stop("No variables in data tables specified by \'targets\' are numeric.")
            else
              coxnetTrainInterface(measurementsTrain, survivalTrain, ...)
          })



################################################################################
#
# Predict Interface
#
################################################################################




# Matrix of numeric measurements.
setGeneric("coxnetPredictInterface", function(model, measurementsTest, ...)
  standardGeneric("coxnetPredictInterface"))

setMethod("coxnetPredictInterface", c("coxnet", "matrix"),
          function(model, measurementsTest, ...)
          {
            coxnetPredictInterface(model, DataFrame(measurementsTest, check.names = FALSE), ...)
          })

# Clinical data only.
setMethod("coxnetPredictInterface", c("coxnet", "DataFrame"), function(model, measurementsTest, survivalTest = NULL, lambda, ..., returnType = c("both", "class", "score"), verbose = 3)
{ # ... just consumes emitted tuning variables from .doTrain which are unused.
  if(!is.null(survivalTest))
  {
    splitDataset <- .splitDataAndOutcomes(measurementsTest, survivalTest)  # Remove any classes, if present.
    measurementsTest <- splitDataset[["measurements"]]
  }
  
  returnType <- match.arg(returnType)
  
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using cox model.")
  
  if(missing(lambda)) # Tuning parameters are not passed to prediction functions.
    lambda <- attr(model, "tune")[["lambda"]] # Sneak it in as an attribute on the model.
  
  testMatrix <- glmnet::makeX(as(measurementsTest, "data.frame"))
  testMatrix <- testMatrix[, rownames(model[["beta"]])]
  
  offset <- attr(model, "tune")[["offset"]]
  model$offset <- TRUE
  
  survPredictions <- predict(model, testMatrix, s = lambda, type = "link", newoffset = offset)
  survScores <- predict(model, testMatrix, s = lambda, type = "response", newoffset = offset)
  
  data.frame(link = survPredictions[, 1], relativeRisk = survScores[, 1], check.names = FALSE)
})

# One or more omics data sets, possibly with clinical data.
setMethod("coxnetPredictInterface", c("coxnet", "MultiAssayExperiment"),
          function(model, measurementsTest, targets = names(measurementsTest), ...)
          {
            tables <- .MAEtoWideTable(measurementsTest, targets)
            measurementsTest <- tables[["dataTable"]]
            
            coxnetPredictInterface(model, measurementsTest, ...)
          })



################################################################################
#
# Get selected features
#
################################################################################
