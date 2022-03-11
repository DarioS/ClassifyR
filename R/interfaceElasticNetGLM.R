################################################################################
#
# Train Interface
#
################################################################################

#' An Interface for glmnet Package's glmnet Function
#' 
#' An elastic net GLM classifier uses a penalty which is a combination of a
#' lasso penalty and a ridge penalty, scaled by a lambda value, to fit a sparse
#' linear model to the data.
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
#' @param lambda The lambda value passed directly to
#' \code{\link[glmnet]{glmnet}} if the training function is used or passed as
#' \code{s} to \code{\link[glmnet]{predict.glmnet}} if the prediction function
#' is used.
#' @param measurementsTest An object of the same class as \code{measurementsTrain} with no
#' samples in common with \code{measurementsTrain} and the same number of features
#' as it.
#' @param targets If \code{measurementsTrain} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"sampleInfo"} is also a valid value
#' and specifies that integer variables from the sample information data table will be
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
#'     genesMatrix <- sapply(1:100, function(sample) rnorm(25, 9, 0.3))
#'     genesMatrix <- rbind(genesMatrix, t(sapply(1:25, function(sample)
#'                                       c(rnorm(75, 9, 0.3), rnorm(25, 14, 0.3)))))
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     rownames(genesMatrix) <- paste("Sample", 1:nrow(genesMatrix))
#'     colnames(genesMatrix) <- paste("Gene", 1:ncol(genesMatrix))
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
#' @rdname elasticNetGLM
#' @export
setGeneric("elasticNetGLMtrainInterface", function(measurementsTrain, ...)
standardGeneric("elasticNetGLMtrainInterface"))

#' @rdname elasticNetGLM
#' @export
setMethod("elasticNetGLMtrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurementsTrain, classesTrain, ...)
{
  elasticNetGLMtrainInterface(S4Vectors::DataFrame(measurementsTrain, check.names = FALSE), classesTrain, ...)
})

# Sample information data or one of the other inputs, transformed.
#' @rdname elasticNetGLM
#' @seealso \code{\link{elasticNetFeatures}} for a function used to extract the features
#' with non-zero coefficients from the model.
#' @export
setMethod("elasticNetGLMtrainInterface", "DataFrame", function(measurementsTrain, classesTrain, lambda = NULL, ..., verbose = 3)
{
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting elastic net regularised GLM classifier to data.")
  
  splitDataset <- .splitDataAndOutcomes(measurementsTrain, classesTrain, restrict = NULL)
  measurementsTrain <- data.frame(splitDataset[["measurements"]], check.names = FALSE)
  measurementsMatrix <- glmnet::makeX(as(measurementsTrain, "data.frame"))

  fitted <- glmnet::glmnet(measurementsMatrix, splitDataset[["outcomes"]], family = "multinomial", ...)

  if(is.null(lambda)) # fitted has numerous models for automatically chosen lambda values.
  { # Pick one lambda based on resubstitution performance.
    bestLambda <- fitted[["lambda"]][which.min(sapply(fitted[["lambda"]], function(lambda) # Largest Lambda with minimum balanced error rate.
    {
      classPredictions <- factor(as.character(predict(fitted, measurementsMatrix, s = lambda, type = "class")), levels = fitted[["classnames"]])
      calcExternalPerformance(splitDataset[["outcomes"]], classPredictions, "Balanced Error")
    }))[1]]
    attr(fitted, "tune") <- list(lambda = bestLambda)
  }
  fitted
})

# One or more omics datasets, possibly with sample information data.
#' @rdname elasticNetGLM
#' @export
setMethod("elasticNetGLMtrainInterface", "MultiAssayExperiment",
function(measurementsTrain, targets = names(measurementsTrain), classesTrain, ...)
{
  tablesAndOutcomes <- .MAEtoWideTable(measurementsTrain, targets, classesTrain, restrict = NULL)
  measurementsTrain <- tablesAndOutcomes[["dataTable"]]
  classesTrain <- tablesAndOutcomes[["outcomes"]]
  
  if(ncol(measurementsTrain) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    elasticNetGLMtrainInterface(measurementsTrain, classesTrain, ...)
})



################################################################################
#
# Predict Interface
#
################################################################################

# Matrix of numeric measurements.
#' @rdname elasticNetGLM
#' @export
setGeneric("elasticNetGLMpredictInterface", function(model, measurementsTest, ...)
standardGeneric("elasticNetGLMpredictInterface"))

#' @rdname elasticNetGLM
#' @export
setMethod("elasticNetGLMpredictInterface", c("multnet", "matrix"),
          function(model, measurementsTest, ...)
{
  elasticNetGLMpredictInterface(model, DataFrame(measurementsTest, check.names = FALSE), ...)
})

# Sample information data, for example.
#' @rdname elasticNetGLM
#' @export
setMethod("elasticNetGLMpredictInterface", c("multnet", "DataFrame"), function(model, measurementsTest, lambda, ..., returnType = c("both", "class", "score"), verbose = 3)
{ # ... just consumes emitted tuning variables from .doTrain which are unused.
  returnType <- match.arg(returnType)
  
  # Ensure that testing data has same columns names in same order as training data.
  # Remove those annoying backquotes which glmnet adds if variables have spaces in names.
  measurementsTest <- measurementsTest[, gsub('`', '', rownames(model[["beta"]][[1]]))]
  
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained elastic net regularised GLM classifier.")

  if(missing(lambda)) # Tuning parameters are not passed to prediction functions.
    lambda <- attr(model, "tune")[["lambda"]] # Sneak it in as an attribute on the model.

  testMatrix <- glmnet::makeX(as(measurementsTest, "data.frame"))
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

# One or more omics data sets, possibly with sample information data.
#' @rdname elasticNetGLM
#' @export
setMethod("elasticNetGLMpredictInterface", c("multnet", "MultiAssayExperiment"),
          function(model, measurementsTest, targets = names(measurementsTest), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurementsTest, targets)
  measurementsTest <- tablesAndClasses[["dataTable"]]

  elasticNetGLMpredictInterface(model, measurementsTest, ...)
})

################################################################################
#
# Get selected features
#
################################################################################

#' Extract Vectors of Ranked and Selected Features From an Elastic Net GLM
#' Object
#' 
#' Provides a ranking of features based on the magnitude of fitted GLM
#' coefficients. Also provides the selected features which are those with a
#' non-zero coefficient.
#' 
#' 
#' @aliases elasticNetFeatures elasticNetFeatures,multnet-method
#' @param model A fitted multinomial GLM which was created by
#' \code{\link[glmnet]{glmnet}}.
#' @return An \code{list} object. The first element is a vector or data frame
#' of ranked features, the second is a vector or data frame of selected
#' features.
#' @author Dario Strbenac
#' @examples
#' 
#'     if(require(glmnet))
#'     {
#'       # Genes 76 to 100 have differential expression.
#'     genesMatrix <- sapply(1:100, function(sample) rnorm(25, 9, 0.3))
#'     genesMatrix <- rbind(genesMatrix, t(sapply(1:25, function(sample)
#'                                       c(rnorm(75, 9, 0.3), rnorm(25, 14, 0.3)))))
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     rownames(genesMatrix) <- paste("Sample", 1:nrow(genesMatrix))
#'     colnames(genesMatrix) <- paste("Gene", 1:ncol(genesMatrix))
#'                                              
#'       # alpha is a user-specified tuning parameter.
#'       # lambda is automatically tuned, based on glmnet defaults, if not user-specified.                 
#'       CVparams <- CrossValParams("k-Fold")
#'       
#'       trainParams <- TrainParams(elasticNetGLMtrainInterface, nlambda = 500)
#'       predictParams <- PredictParams(elasticNetGLMpredictInterface)
#'       modParams <- ModellingParams(selectParams = NULL, trainParams = trainParams,
#'                                    predictParams = predictParams)
#'                                    
#'       classified <- runTests(genesMatrix, classes, CVparams, modParams)
#'                                         
#'       elasticNetFeatures(models(classified)[[1]])
#'     }
#' @rdname elasticNetFeatures
#' @export
#' @usage NULL
setGeneric("elasticNetFeatures", function(model, ...)
  standardGeneric("elasticNetFeatures"))

#' @rdname elasticNetFeatures
#' @export
setMethod("elasticNetFeatures", "multnet",
          function(model)
          {
            inputFeatures <- rownames(model[["beta"]][[1]])            
            # Floating point numbers test for equality.
            whichCoefficientColumn <- which(abs(model[["lambda"]] - attr(model, "tune")[["lambda"]]) < 0.00001)[1]
            coefficientsUsed <- sapply(model[["beta"]], function(classCoefficients) classCoefficients[, whichCoefficientColumn])
            featureScores <- rowSums(abs(coefficientsUsed))
            rankedFeatures <- inputFeatures[order(featureScores, decreasing = TRUE)]
            selectedFeatures <- inputFeatures[featureScores != 0]
            
            # Colon is a reserved symbol for separating data name and feature name, which is necessary for identifiability of MultiAssayExperiment features. It is not permitted in feature names.
            if(grepl(':', selectedFeatures[1]) == TRUE) # Convert to data.frame.
            {
              selectedFeatures <- do.call(rbind, strsplit(selectedFeatures, ':'))
              rankedFeatures <- do.call(rbind, strsplit(rankedFeatures, ':'))
              colnames(selectedFeatures) <- colnames(rankedFeatures) <- c("dataset", "feature")
            }
            list(unique(rankedFeatures), selectedFeatures)
          })