################################################################################
#
# Train Interface
#
################################################################################

#' An Interface for Ordinary Generalised Linear Models (GLMs) with Binary Outcome
#' 
#' No feature selection is done, so some should be done beforehand or otherwise the
#' elastic net GLM classifier should be used instead of this classifier.
#' 
#' The value of the \code{family} parameter is fixed to \code{binomal}.
#' 
#' @name Ordinary GLM Interface
#' @aliases GLMinterface GLMtrainInterface
#' GLMpredictInterface GLMtrainInterface,matrix-method
#' GLMtrainInterface,DataFrame-method
#' GLMtrainInterface,MultiAssayExperiment-method
#' GLMpredictInterface,glm,matrix-method
#' GLMpredictInterface,glm,DataFrame-method
#' GLMpredictInterface,glm,MultiAssayExperiment-method
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
#' options that are used by the \code{\link{glm}} function. For the testing
#' function, this variable simply contains any parameters passed from the
#' classification framework to it which aren't used by the \code{\link{predict.glm}}
#' function.
#' @param model A trained GLM, as created by the \code{glm} function.
#' @param returnType Default: \code{"both"}. Either \code{"class"},
#' \code{"score"} or \code{"both"}.  Sets the return value from the prediction
#' to either a vector of class labels, matrix of scores for each class, or both
#' labels and scores in a \code{data.frame}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return For \code{GLMtrainInterface}, an object of type
#' \code{list}. For \code{GLMpredictInterface}, either a factor
#' vector of predicted classes, a matrix of scores for each class, or a table
#' of both the class labels and class scores, depending on the setting of
#' \code{returnType}.
#' @author Dario Strbenac
#' @examples
#' 
#'     # Genes 76 to 100 have differential expression.
#'     genesMatrix <- sapply(1:100, function(sample) rnorm(25, 9, 0.3))
#'     genesMatrix <- rbind(genesMatrix, t(sapply(1:25, function(sample)
#'                                       c(rnorm(75, 9, 0.3), rnorm(25, 14, 0.3)))))
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     rownames(genesMatrix) <- paste("Sample", 1:nrow(genesMatrix))
#'     colnames(genesMatrix) <- paste("Gene", 1:ncol(genesMatrix))
#'     
#'     CVparams <- CrossValParams("k-Fold")
#'     trainParams <- TrainParams(GLMtrainInterface)
#'     predictParams <- PredictParams(GLMpredictInterface)
#'     modParams <- ModellingParams(selectParams = NULL, trainParams = trainParams,
#'                                    predictParams = predictParams)
#'     classified <- runTests(genesMatrix[, 76:80], classes, CVparams, modParams)
#'                            
#'     classified <- calcCVperformance(classified, "Balanced Error")
#'     head(tunedParameters(classified))
#'     performance(classified)
#' @usage NULL
#' @rdname GLM
#' @export
setGeneric("GLMtrainInterface", function(measurementsTrain, ...)
standardGeneric("GLMtrainInterface"))

#' @rdname GLM
#' @export
setMethod("GLMtrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurementsTrain, classesTrain, ...)
{
  GLMtrainInterface(S4Vectors::DataFrame(measurementsTrain, check.names = FALSE), classesTrain, ...)
})

# Sample information data or one of the other inputs, transformed.
#' @rdname GLM
#' @export
setMethod("GLMtrainInterface", "DataFrame", function(measurementsTrain, classesTrain, ..., verbose = 3)
{
  if(verbose == 3)
    message("Fitting GLM classifier to data.")
    
    if(is.factor(classesTrain))
    {
        fitData <- cbind(measurementsTrain, class = classesTrain)
        classesTrain <- "class" # Column name for glm fit.
    } else {fitData <- measurementsTrain}
  glm(class ~ . + 0, family = binomial, data = fitData, ...)
})

# One or more omics assays, possibly with sample information data.
#' @rdname GLM
#' @export
setMethod("GLMtrainInterface", "MultiAssayExperiment",
function(measurementsTrain, targets = names(measurementsTrain), classesTrain, ...)
{
  tablesAndOutcomes <- .MAEtoWideTable(measurementsTrain, targets, classesTrain, restrict = NULL)
  measurementsTrain <- tablesAndOutcomes[["dataTable"]]
  classesTrain <- tablesAndOutcomes[["outcomes"]]
  
  if(ncol(measurementsTrain) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    GLMtrainInterface(measurementsTrain, classesTrain, ...)
})



################################################################################
#
# Predict Interface
#
################################################################################

# Matrix of numeric measurements.
#' @usage NULL
#' @export
setGeneric("GLMpredictInterface", function(model, measurementsTest, ...)
standardGeneric("GLMpredictInterface"))

#' @rdname GLM
#' @export
setMethod("GLMpredictInterface", c("glm", "matrix"),
          function(model, measurementsTest, ...)
{
  GLMpredictInterface(model, DataFrame(measurementsTest, check.names = FALSE), ...)
})

#' @rdname GLM
#' @export
setMethod("GLMpredictInterface", c("glm", "DataFrame"), function(model, measurementsTest, returnType = c("both", "class", "score"), verbose = 3)
{
  returnType <- match.arg(returnType)
  
  if(verbose == 3)
    message("Predicting classes using trained GLM classifier.")
  
  predictions <- predict(model, measurementsTest, type = "response")
  classes <- levels(model[["model"]][["class"]])
  classPredictions <- factor(ifelse(predictions >= 0.5, classes[2], classes[1]), classes)
  classScores <- matrix(c(1 - predictions, predictions), ncol = 2)
  colnames(classScores) <- classes
  
  switch(returnType, class = classPredictions, # Factor vector.
         score = classScores, # Numeric matrix.
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
})

# One or more omics data sets, possibly with sample information data.
#' @rdname GLM
#' @export
setMethod("GLMpredictInterface", c("glm", "MultiAssayExperiment"),
          function(model, measurementsTest, targets = names(measurementsTest), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurementsTest, targets)
  measurementsTest <- tablesAndClasses[["dataTable"]]

  GLMpredictInterface(model, measurementsTest, ...)
})