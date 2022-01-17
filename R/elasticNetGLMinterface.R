#' An Interface for glmnet Package's glmnet Function
#' 
#' An elastic net GLM classifier uses a penalty which is a combination of a
#' lasso penalty and a ridge penalty, scaled by a lambda value, to fit a sparse
#' linear model to the data.
#' 
#' If \code{measurements} is an object of class \code{MultiAssayExperiment},
#' the factor of sample classes must be stored in the DataFrame accessible by
#' the \code{colData} function with column name \code{"class"}.
#' 
#' The value of the \code{family} parameter is fixed to \code{"multinomial"} so
#' that classification with more than 2 classes is possible and
#' \code{type.multinomial} is fixed to \code{"grouped"} so that a grouped lasso
#' penalty is used. During classifier training, if more than one lambda value
#' is considered by specifying a vector of them as input or leaving the default
#' value of NULL, then the chosen value is determined based on classifier
#' resubstitution error rate.
#' 
#' @aliases elasticNetGLMinterface elasticNetGLMtrainInterface
#' elasticNetGLMpredictInterface elasticNetGLMtrainInterface,matrix-method
#' elasticNetGLMtrainInterface,DataFrame-method
#' elasticNetGLMtrainInterface,MultiAssayExperiment-method
#' elasticNetGLMpredictInterface,multnet,matrix-method
#' elasticNetGLMpredictInterface,multnet,DataFrame-method
#' elasticNetGLMpredictInterface,multnet,MultiAssayExperiment-method
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
#' @param model A trained elastic net GLM, as created by the \code{glmnet}
#' function.
#' @param returnType Default: \code{"both"}. Either \code{"class"},
#' \code{"score"} or \code{"both"}.  Sets the return value from the prediction
#' to either a vector of class labels, matrix of scores for each class, or both
#' labels and scores in a \code{data.frame}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return For \code{elasticNetGLMtrainInterface}, an object of type
#' \code{glmnet}. For \code{elasticNetGLMpredictInterface}, either a factor
#' vector of predicted classes, a matrix of scores for each class, or a table
#' of both the class labels and class scores, depending on the setting of
#' \code{returnType}.
#' @author Dario Strbenac
#' @examples
#' 
#'   if(require(glmnet))
#'   {
#'     # Genes 76 to 100 have differential expression.
#'     genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
#'     genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
#'                                       c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix))
#'     rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
#'     
#'     CVparams <- CrossValParams("k-Fold")
#'       
#'     trainParams <- TrainParams(elasticNetGLMtrainInterface, nlambda = 500)
#'     predictParams <- PredictParams(elasticNetGLMpredictInterface)
#'     modParams <- ModellingParams(selectParams = NULL, trainParams = trainParams,
#'                                    predictParams = predictParams)
#'     classified <- runTests(genesMatrix, classes, CVparams, modParams)
#'                            
#'     classified <- calcCVperformance(classified, "Balanced Error")
#'     head(tunedParameters(classified))
#'     performance(classified)
#'   }
#'   
#' @export
setGeneric("elasticNetGLMtrainInterface", function(measurements, ...)
standardGeneric("elasticNetGLMtrainInterface"))

setMethod("elasticNetGLMtrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  elasticNetGLMtrainInterface(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

# Clinical data or one of the other inputs, transformed.
setMethod("elasticNetGLMtrainInterface", "DataFrame", function(measurements, classes, lambda = NULL, ..., verbose = 3)
{
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting elastic net regularised GLM classifier to data.")
  
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- data.frame(splitDataset[["measurements"]], check.names = FALSE)
  measurementsMatrix <- glmnet::makeX(as(measurements, "data.frame"))

  fitted <- glmnet::glmnet(measurementsMatrix, splitDataset[["classes"]], family = "multinomial", ...)

  if(is.null(lambda)) # fitted has numerous models for automatically chosen lambda values.
  { # Pick one lambda based on resubstitution performance.
    bestLambda <- fitted[["lambda"]][which.min(sapply(fitted[["lambda"]], function(lambda) # Largest Lambda with minimum balanced error rate.
    {
      classPredictions <- factor(as.character(predict(fitted, measurementsMatrix, s = lambda, type = "class")), levels = fitted[["classnames"]])
      calcExternalPerformance(classes, classPredictions, "Balanced Error")
    }))[1]]
    attr(fitted, "tune") <- list(lambda = bestLambda)
  }
  fitted
})

# One or more omics datasets, possibly with clinical data.
setMethod("elasticNetGLMtrainInterface", "MultiAssayExperiment",
function(measurements, targets = names(measurements), classes, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  
  if(ncol(measurements) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    elasticNetGLMtrainInterface(measurements, classes, ...)
})

# Matrix of numeric measurements.
setGeneric("elasticNetGLMpredictInterface", function(model, test, ...)
standardGeneric("elasticNetGLMpredictInterface"))

setMethod("elasticNetGLMpredictInterface", c("multnet", "matrix"),
          function(model, test, ...)
{
  elasticNetGLMpredictInterface(model, DataFrame(t(test), check.names = FALSE), ...)
})

# Clinical data only.
setMethod("elasticNetGLMpredictInterface", c("multnet", "DataFrame"), function(model, test, classes = NULL, lambda, ..., returnType = c("both", "class", "score"), verbose = 3)
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
    message("Predicting classes using trained elastic net regularised GLM classifier.")

  if(missing(lambda)) # Tuning parameters are not passed to prediction functions.
    lambda <- attr(model, "tune")[["lambda"]] # Sneak it in as an attribute on the model.

  testMatrix <- glmnet::makeX(as(test, "data.frame"))
  testMatrix <- testMatrix[, rownames(model[["beta"]][[1]])]
  
  classPredictions <- factor(as.character(predict(model, testMatrix, s = lambda, type = "class")), levels = model[["classnames"]])
  classScores <- predict(model, testMatrix, s = lambda, type = "response")[, , 1]
  
  if(is.matrix(classScores))
    classScores <- classScores[, model[["classnames"]]]
  else # Leave-one-out cross-validation likely used and glmnet doesn't have consistent return types.
    classScores <- t(classScores[model[["classnames"]]])
  
  switch(returnType, class = classPredictions, # Factor vector.
         score = classScores, # Numeric matrix.
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
})

# One or more omics data sets, possibly with clinical data.
setMethod("elasticNetGLMpredictInterface", c("multnet", "MultiAssayExperiment"),
          function(model, test, targets = names(test), ...)
{
  tablesAndClasses <- .MAEtoWideTable(test, targets)
  test <- tablesAndClasses[["dataTable"]]

  elasticNetGLMpredictInterface(model, test, ...)
})
