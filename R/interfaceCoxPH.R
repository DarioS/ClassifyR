################################################################################
#
# Train interface
#
################################################################################


#' An Interface for survival Package's coxph Function
#' 
#' 
#' If \code{measurements} is an object of class \code{MultiAssayExperiment},
#' the factor of sample classes must be stored in the DataFrame accessible by
#' the \code{colData} function with column name \code{"class"}.
#' 
#' @aliases coxphInterface coxphTrainInterface
#' coxphInterface,matrix-method
#' coxphInterface,DataFrame-method
#' coxphInterface,MultiAssayExperiment-method
#' coxphPredictInterface
#' coxphPredictInterface,coxph,matrix-method
#' coxphPredictInterface,coxph,DataFrame-method
#' coxphPredictInterface,coxph,MultiAssayExperiment-method
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix}, the rows are features, and the columns are samples.
#' @param classes Either a vector of class labels of class \code{\link{factor}}
#' of the same length as the number of samples in \code{measurements} or if the
#' measurements are of class \code{DataFrame} a character vector of length 1
#' containing the column name in \code{measurement} is also permitted. Not used
#' if \code{measurements} is a \code{MultiAssayExperiment} object.
#' @param model A trained random coxph classifier, as created by
#' \code{coxphTrainInterface}, which has the same form as the output of
#' \code{\link[coxph]{coxph}}.
#' @param test An object of the same class as \code{measurements} with no
#' samples in common with \code{measurements} and the same number of features
#' as it.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"clinical"} is also a valid value
#' and specifies that integer variables from the clinical data table will be
#' used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method (e.g. \code{verbose}) or options which are accepted
#' by the \code{\link[coxph]{coxph}} or
#' \code{\link[coxph]{predict.coxph}} functions.
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
#' @importFrom survival coxph predict concordance
#' @export
setGeneric("coxphTrainInterface", function(measurements, ...)
standardGeneric("coxphTrainInterface"))

setMethod("coxphTrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  coxphTrainInterface(DataFrame(t(measurements), check.names = FALSE),
                             classes, ...)
})

# Clinical data or one of the other inputs, transformed.
setMethod("coxphTrainInterface", "DataFrame", function(measurements, classes, ..., verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)

  if(!requireNamespace("survival", quietly = TRUE))
    stop("The package 'survival' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting coxph classifier to training data and making predictions on test
            data.")
  
  classes <- splitDataset[["classes"]]
  measurements <- splitDataset[["measurements"]]
  
  survival::coxph(classes~., measurements)
})

setMethod("coxphTrainInterface", "MultiAssayExperiment",
function(measurements, targets = names(measurements), classes, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes, restrict = NULL)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  
  coxphTrainInterface(measurements, classes, ...)
})



################################################################################
#
# Predict Interface
#
################################################################################



#' @export
setGeneric("coxphPredictInterface", function(model, test, ...)
  standardGeneric("coxphPredictInterface"))

setMethod("coxphPredictInterface", c("coxph", "matrix"), function(model, test, ...)
{
  coxphPredictInterface(forest, DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("coxphPredictInterface", c("coxph", "DataFrame"),
function(model, test, ..., verbose = 3)
{
  
  output <- sapply(c("lp", "risk"), function(type)predict(model, as.data.frame(test), type = type), simplify = TRUE)
  output <- data.frame(output)
  output$class <- output[, "lp"]
  output

})

# One or more omics data sets, possibly with clinical data.
setMethod("coxphPredictInterface", c("coxph", "MultiAssayExperiment"),
          function(model, test, targets = names(test), ...)
{
  testingTable <- .MAEtoWideTable(test, targets)
  coxphPredictInterface(model, testingTable, ...)
})


