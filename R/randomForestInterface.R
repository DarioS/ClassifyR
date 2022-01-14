#' An Interface for randomForest Package's randomForest Function
#' 
#' A random forest classifier builds multiple decision trees and uses the
#' predictions of the trees to determine a single prediction for each test
#' sample.
#' 
#' If \code{measurements} is an object of class \code{MultiAssayExperiment},
#' the factor of sample classes must be stored in the DataFrame accessible by
#' the \code{colData} function with column name \code{"class"}.
#' 
#' @aliases randomForestInterface randomForestTrainInterface
#' randomForestTrainInterface,matrix-method
#' randomForestTrainInterface,DataFrame-method
#' randomForestTrainInterface,MultiAssayExperiment-method
#' randomForestPredictInterface
#' randomForestPredictInterface,randomForest,matrix-method
#' randomForestPredictInterface,randomForest,DataFrame-method
#' randomForestPredictInterface,randomForest,MultiAssayExperiment-method
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix}, the rows are features, and the columns are samples.
#' @param classes Either a vector of class labels of class \code{\link{factor}}
#' of the same length as the number of samples in \code{measurements} or if the
#' measurements are of class \code{DataFrame} a character vector of length 1
#' containing the column name in \code{measurement} is also permitted. Not used
#' if \code{measurements} is a \code{MultiAssayExperiment} object.
#' @param forest A trained random forest classifier, as created by
#' \code{randomForestTrainInterface}, which has the same form as the output of
#' \code{\link[randomForest]{randomForest}}.
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
#' by the \code{\link[randomForest]{randomForest}} or
#' \code{\link[randomForest]{predict.randomForest}} functions.
#' @param returnType Default: \code{"both"}. Either \code{"class"},
#' \code{"score"} or \code{"both"}.  Sets the return value from the prediction
#' to either a vector of class labels, score for a sample belonging to the
#' second class, as determined by the factor levels, or both labels and scores
#' in a \code{data.frame}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return For \code{randomForestTrainInterface}, the trained random forest.
#' For \code{randomForestPredictInterface}, either a factor vector of predicted
#' classes, a matrix of scores for each class, or a table of both the class
#' labels and class scores, depending on the setting of \code{returnType}.
#' @author Dario Strbenac
#' @examples
#' 
#'   if(require(randomForest))
#'   {
#'     # Genes 76 to 100 have differential expression.
#'     genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
#'     genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
#'                                       c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix), sep = ' ')
#'     rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix), sep = '-')
#'     trainingSamples <- c(1:20, 26:45)
#'     testingSamples <- c(21:25, 46:50)
#'     
#'     trained <- randomForestTrainInterface(genesMatrix[, trainingSamples],
#'                                           classes[trainingSamples])
#'     predicted <- randomForestPredictInterface(trained, genesMatrix[, testingSamples])
#'   }
#' 
#' @export
setGeneric("randomForestTrainInterface", function(measurements, ...)
standardGeneric("randomForestTrainInterface"))

#' @export
setGeneric("randomForestPredictInterface", function(models, test, ...)
           standardGeneric("randomForestPredictInterface"))

setMethod("randomForestTrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  randomForestTrainInterface(DataFrame(t(measurements), check.names = FALSE),
                             classes, ...)
})

# Clinical data or one of the other inputs, transformed.
setMethod("randomForestTrainInterface", "DataFrame", function(measurements, classes, ..., verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)

  if(!requireNamespace("randomForest", quietly = TRUE))
    stop("The package 'randomForest' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting random forest classifier to training data and making predictions on test
            data.")

  # Convert to base data.frame as randomForest doesn't understand DataFrame.
  randomForest::randomForest(as(splitDataset[["measurements"]], "data.frame"), splitDataset[["classes"]], keep.forest = TRUE, ...)
})

setMethod("randomForestTrainInterface", "MultiAssayExperiment",
function(measurements, targets = names(measurements), classes, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes, restrict = NULL)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  
  randomForestTrainInterface(measurements, classes, ...)
})

setGeneric("randomForestPredictInterface", function(forest, test, ...)
           standardGeneric("randomForestPredictInterface"))

setMethod("randomForestPredictInterface", c("randomForest", "matrix"), function(forest, test, ...)
{
  randomForestPredictInterface(forest, DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("randomForestPredictInterface", c("randomForest", "DataFrame"),
function(forest, test, ..., returnType = c("both", "class", "score"), verbose = 3)
{
  returnType <- match.arg(returnType)
  if(verbose == 3)
    message("Predicting using random forest.")  
  
  classPredictions <- predict(forest, test)
  classScores <- predict(forest, test, type = "vote")[, forest[["classes"]], drop = FALSE]
  switch(returnType, class = classPredictions,
         score = classScores,
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
})

# One or more omics data sets, possibly with clinical data.
setMethod("randomForestPredictInterface", c("randomForest", "MultiAssayExperiment"),
          function(forest, test, targets = names(test), ...)
{
  testingTable <- .MAEtoWideTable(test, targets)
  randomForestPredictInterface(models, testingTable, ...)
})