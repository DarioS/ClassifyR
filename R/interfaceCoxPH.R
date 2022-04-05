################################################################################
#
# Train interface
#
################################################################################


#' An Interface for survival Package's coxph Function
#' 
#' Cox proportional hazards.
#' 
#' @aliases coxphInterface coxphTrainInterface
#' coxphInterface,matrix-method
#' coxphInterface,DataFrame-method
#' coxphInterface,MultiAssayExperiment-method
#' coxphPredictInterface
#' coxphPredictInterface,coxph,matrix-method
#' coxphPredictInterface,coxph,DataFrame-method
#' coxphPredictInterface,coxph,MultiAssayExperiment-method
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
#' @param model A trained coxph classifier, as created by
#' \code{coxphTrainInterface}, which has the same form as the output of
#' \code{\link[survival]{coxph}}.
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
#' by the \code{\link[survival]{coxph}} or \code{\link[survival]{predict}} functions.
#' @param returnType Default: \code{"both"}. Either \code{"class"},
#' \code{"score"} or \code{"both"}.  Sets the return value from the prediction
#' to either a vector of class labels, score for a sample belonging to the
#' second class, as determined by the factor levels, or both labels and scores
#' in a \code{data.frame}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return For \code{coxphTrainInterface}, the trained coxph.
#' For \code{coxphPredictInterface}, either a factor vector of predicted
#' classes, a matrix of scores for each class, or a table of both the class
#' labels and class scores, depending on the setting of \code{returnType}.
#' @examples
#' #' 
#'   # if(require(randomForest))
#'   # {
#'   #   # Genes 76 to 100 have differential expression.
#'   #   genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
#'   #   genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
#'   #                                     c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
#'   #   classes <- factor(rep(c("Poor", "Good"), each = 25))
#'   #   colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix), sep = ' ')
#'   #   rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix), sep = '-')
#'   #   trainingSamples <- c(1:20, 26:45)
#'   #   testingSamples <- c(21:25, 46:50)
#'   # 
#'   #   trained <- randomForestTrainInterface(genesMatrix[, trainingSamples],
#'   #                                         classes[trainingSamples])
#'   #   predicted <- randomForestPredictInterface(trained, genesMatrix[, testingSamples])
#'   # }
#' 
#' @importFrom survival coxph concordance
#' @export
setGeneric("coxphTrainInterface", function(measurementsTrain, ...) standardGeneric("coxphTrainInterface"))

setMethod("coxphTrainInterface", "matrix", function(measurementsTrain, survivalTrain, ...)
{
  coxphTrainInterface(S4Vectors::DataFrame(measurementsTrain, check.names = FALSE), survivalTrain, ...)
})

# Clinical data or one of the other inputs, transformed.
setMethod("coxphTrainInterface", "DataFrame", function(measurementsTrain, survivalTrain, ..., verbose = 3)
{
  if(!requireNamespace("survival", quietly = TRUE))
    stop("The package 'survival' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting coxph classifier to training data and making predictions on test
            data.")

  splitDataset <- .splitDataAndOutcomes(measurementsTrain, survivalTrain)  
  survivalTrain <- splitDataset[["outcomes"]]
  measurementsTrain <- splitDataset[["measurements"]]
  
  survival::coxph(survivalTrain ~ ., measurementsTrain)
})

setMethod("coxphTrainInterface", "MultiAssayExperiment", function(measurementsTrain, targets = names(measurementsTrain), survivalTrain, ...)
{
  tablesAndSurvival <- .MAEtoWideTable(measurementsTrain, targets, survivalTrain, restrict = NULL)
  measurementsTrain <- tablesAndSurvival[["dataTable"]]
  survivalTrain <- tablesAndSurvival[["outcomes"]]
  
  coxphTrainInterface(measurementsTrain, survivalTrain, ...)
})



################################################################################
#
# Predict Interface
#
################################################################################



#' @export
setGeneric("coxphPredictInterface", function(model, measurementsTest, ...)
  standardGeneric("coxphPredictInterface"))

setMethod("coxphPredictInterface", c("coxph", "matrix"), # Matrix of numeric measurements.
          function(model, measurementsTest, ...)
{
  coxphPredictInterface(model, S4Vectors::DataFrame(measurementsTest, check.names = FALSE), ...)
})

setMethod("coxphPredictInterface", c("coxph", "DataFrame"),
function(model, measurementsTest, ..., verbose = 3)
{
  # Linear predictor ("lp") is default prediction format.
  predict(model, as.data.frame(measurementsTest))
})

# One or more omics data sets, possibly with clinical data.
setMethod("coxphPredictInterface", c("coxph", "MultiAssayExperiment"),
          function(model, measurementsTest, targets = names(measurementsTest), ...)
{
  testingTable <- .MAEtoWideTable(measurementsTest, targets)
  coxphPredictInterface(model, testingTable, ...)
})


