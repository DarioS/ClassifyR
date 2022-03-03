#' An Interface for e1071 Package's Support Vector Machine Classifier.
#' 
#' \code{SVMtrainInterface} generates a trained SVM classifier and
#' \code{SVMpredictInterface} uses it to make predictions on a test data set.
#' 
#' @aliases SVMtrainInterface SVMtrainInterface,matrix-method
#' SVMtrainInterface,DataFrame-method
#' SVMtrainInterface,MultiAssayExperiment-method SVMpredictInterface
#' SVMpredictInterface,svm,matrix-method
#' SVMpredictInterface,svm,DataFrame-method
#' SVMpredictInterface,svm,MultiAssayExperiment-method

#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
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
#' @param measurementsTest An object of the same class as \code{measurementsTrain} with no
#' samples in common with \code{measurementsTrain} and the same number of features
#' as it.
#' @param returnType Default: \code{"both"}. Either \code{"class"},
#' \code{"score"} or \code{"both"}. Sets the return value from the prediction
#' to either a vector of class labels, score for a sample belonging to the
#' second class, as determined by the factor levels, or both labels and scores
#' in a \code{data.frame}.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"sampleInfo"} is also a valid value
#' and specifies that numeric variables from the sample information data table will be
#' used.
#' @param model A fitted model as returned by \code{SVMtrainInterface}.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method (e.g. \code{verbose}) or options that are used by
#' the \code{\link[e1071]{svm}} function.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return For \code{SVMtrainInterface}, a trained SVM classifier of type
#' \code{svm}.  For \code{SVMpredictInterface}, either a factor vector of
#' predicted classes, a vector of secores for the second class, or a table of
#' both the class labels and second class scores, depending on the setting of
#' \code{returnType}.
#' @author Dario Strbenac
#' @examples
#' 
#'   if(require(e1071))
#'   {
#'     # Genes 76 to 100 have differential expression.
#'     genesMatrix <- sapply(1:100, function(sample) rnorm(25, 9, 0.3))
#'     genesMatrix <- rbind(genesMatrix, t(sapply(1:25, function(sample)
#'                                       c(rnorm(75, 9, 0.3), rnorm(25, 14, 0.3)))))
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     rownames(genesMatrix) <- paste("Sample", 1:nrow(genesMatrix))
#'     colnames(genesMatrix) <- paste("Gene", 1:ncol(genesMatrix))
#'     trainingSamples <- c(1:20, 26:45)
#'     testingSamples <- c(21:25, 46:50)
#'     
#'     classifier <- SVMtrainInterface(genesMatrix[trainingSamples, ],
#'                                      classes[trainingSamples], kernel = "linear")
#'     SVMpredictInterface(classifier, genesMatrix[testingSamples, ])
#'   }
#' 
#' @export
setGeneric("SVMtrainInterface", function(measurementsTrain, ...)
standardGeneric("SVMtrainInterface"))

setMethod("SVMtrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurementsTrain, classesTrain, ...)
{
  SVMtrainInterface(DataFrame(measurementsTrain, check.names = FALSE), classesTrain, ...)
})

# Sample information data or one of the other inputs, transformed.
setMethod("SVMtrainInterface", "DataFrame", function(measurementsTrain, classesTrain, ..., verbose = 3)
{
  splitDataset <- .splitDataAndOutcomes(measurementsTrain, classesTrain)
  # Classifier requires matrix input data type.
  trainingMatrix <- as.matrix(splitDataset[["measurements"]])

  if(!requireNamespace("e1071", quietly = TRUE))
    stop("The package 'e1071' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting SVM classifier to data.")

  e1071::svm(trainingMatrix, classesTrain, probability = TRUE, ...)
})

setMethod("SVMtrainInterface", "MultiAssayExperiment",
function(measurementsTrain, targets = names(measurementsTrain), classesTrain, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurementsTrain, targets, classesTrain)
  measurementsTrain <- tablesAndClasses[["dataTable"]]
  classesTrain <- tablesAndClasses[["outcomes"]]
  
  if(ncol(measurementsTrain) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    SVMtrainInterface(measurementsTrain, classesTrain, ...)
})


#' @export
setGeneric("SVMpredictInterface", function(model, measurementsTest, ...)
standardGeneric("SVMpredictInterface"))

setMethod("SVMpredictInterface", c("svm", "matrix"),
          function(model, measurementsTest, ...)
{
  SVMpredictInterface(model, DataFrame(measurementsTest, check.names = FALSE), ...)
})

setMethod("SVMpredictInterface", c("svm", "DataFrame"), function(model, measurementsTest, returnType = c("both", "class", "score"), verbose = 3)
{
  returnType <- match.arg(returnType)

  if(!requireNamespace("e1071", quietly = TRUE))
    stop("The package 'e1071' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained SVM classifier.")
  
  # Prediction function depends on test data having same set of columns in same order as
  # selected features used for training.
  colnames(measurementsTest) <- make.names(colnames(measurementsTest)) # svm silently converts feature names.
  measurementsTest <- measurementsTest[, colnames(model[["SV"]])]
  classPredictions <- predict(model, measurementsTest, probability = TRUE)
  
  # e1071 uses attributes to pass back probabilities. Make them a standalone variable.
  classScores <- attr(classPredictions, "probabilities")[, model[["levels"]], drop = FALSE]
  attr(classPredictions, "probabilities") <- NULL
  
  switch(returnType, class = classPredictions, score = classScores,
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
})

setMethod("SVMpredictInterface", c("svm", "MultiAssayExperiment"),
          function(model, measurementsTest, targets = names(measurementsTest), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurementsTest, targets)
  measurementsTest <- tablesAndClasses[["dataTable"]] # Remove any classes, if present.
            
  if(ncol(measurementsTest) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    SVMpredictInterface(model, measurementsTest, ...)
})