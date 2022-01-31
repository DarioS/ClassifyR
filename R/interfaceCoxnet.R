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
#' If \code{measurements} is an object of class \code{MultiAssayExperiment},
#' the factor of sample classes must be stored in the DataFrame accessible by
#' the \code{colData} function with column name \code{"class"}.
#' 
#' The value of the \code{family} parameter is fixed to \code{"cox"} so
#' that classification with survival is possible and
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
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix}, the rows are features, and the columns are samples.  If of
#' type \code{\link{DataFrame}}, the data set is subset to only those features
#' of type \code{integer}.
#' @param classes A vector of class labels of class \code{\link{factor}} of the
#' same length as the number of samples in \code{measurements} if it is a
#' \code{\link{matrix}} (i.e. number of columns) or a \code{\link{DataFrame}}
#' (i.e. number of rows) or a character vector of length 1 containing the
#' column name in \code{measurements} if it is a \code{\link{DataFrame}} or the
#' column name in \code{colData(measurements)} if \code{measurements} is a
#' \code{\link{MultiAssayExperiment}}. If a column name, that column will be
#' removed before training.
#' @param lambda The lambda value passed directly to
#' \code{\link[glmnet]{glmnet}} if the training function is used or passed as
#' \code{s} to \code{\link[glmnet]{predict.glmnet}} if the prediction function
#' is used.
#' @param test An object of the same class as \code{measurements} with no
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
setGeneric("coxnetTrainInterface", function(measurements, ...)
standardGeneric("coxnetTrainInterface"))

setMethod("coxnetTrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  coxnetTrainInterface(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

# Clinical data or one of the other inputs, transformed.
setMethod("coxnetTrainInterface", "DataFrame", function(measurements, classes, lambda = NULL, ..., verbose = 3)
{
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting coxnet model to data.")
  
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- data.frame(splitDataset[["measurements"]], check.names = FALSE)
  measurementsMatrix <- glmnet::makeX(as(measurements, "data.frame"))

  fit <- glmnet::cv.glmnet(measurementsMatrix, splitDataset[["classes"]], family = "cox", type = "C", ...)
  fitted <- fit$glmnet.fit
  
  offset <- -mean(predict(fitted, measurementsMatrix, s = fit$lambda.min, type = "link"))

  attr(fitted, "tune") <- list(lambda = fit$lambda.min, offset = offset)
  fitted
})

# One or more omics datasets, possibly with clinical data.
setMethod("coxnetTrainInterface", "MultiAssayExperiment",
function(measurements, targets = names(measurements), classes, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  
  if(ncol(measurements) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    coxnetTrainInterface(measurements, classes, ...)
})



################################################################################
#
# Predict Interface
#
################################################################################




# Matrix of numeric measurements.
setGeneric("coxnetPredictInterface", function(model, test, ...)
standardGeneric("coxnetPredictInterface"))

setMethod("coxnetPredictInterface", c("coxnet", "matrix"),
          function(model, test, ...)
{
  coxnetPredictInterface(model, DataFrame(t(test), check.names = FALSE), ...)
})

# Clinical data only.
setMethod("coxnetPredictInterface", c("coxnet", "DataFrame"), function(model, test, classes = NULL, lambda, ..., returnType = c("both", "class", "score"), verbose = 3)
{ # ... just consumes emitted tuning variables from .doTrain which are unused.
  if(!is.null(classes))
  {
    splitDataset <- .splitDataAndClasses(test, classes)  # Remove any classes, if present.
    test <- splitDataset[["measurements"]]
  }
  
  returnType <- match.arg(returnType)
  
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using cox model.")

  if(missing(lambda)) # Tuning parameters are not passed to prediction functions.
    lambda <- attr(model, "tune")[["lambda"]] # Sneak it in as an attribute on the model.

  testMatrix <- glmnet::makeX(as(test, "data.frame"))
  testMatrix <- testMatrix[, rownames(model[["beta"]])]
  
  offset <- attr(model, "tune")[["offset"]]
  model$offset <- TRUE
  
  classPredictions <- predict(model, testMatrix, s = lambda, type = "link", newoffset = offset)
  classScores <- predict(model, testMatrix, s = lambda, type = "response", newoffset = offset)
  
  data.frame(link = classPredictions[,1], relativeRisk = classScores[,1], class = classPredictions[,1], check.names = FALSE)
})

# One or more omics data sets, possibly with clinical data.
setMethod("coxnetPredictInterface", c("coxnet", "MultiAssayExperiment"),
          function(model, test, targets = names(test), ...)
{
  tablesAndClasses <- .MAEtoWideTable(test, targets)
  test <- tablesAndClasses[["dataTable"]]

  coxnetPredictInterface(model, test, ...)
})



################################################################################
#
# Get selected features
#
################################################################################



