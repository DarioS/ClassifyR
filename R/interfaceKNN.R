#' An Interface for class Package's knn Function
#' 
#' More details of k Nearest Neighbours are available in the documentation of
#' \code{\link[class]{knn}}.
#' 
#' Data tables which consist entirely of non-numeric data cannot be analysed.
#' If \code{measurements} is an object of class \code{MultiAssayExperiment},
#' the factor of sample classes must be stored in the DataFrame accessible by
#' the \code{colData} function with column name \code{"class"}.
#' 
#' @aliases kNNinterface kNNinterface,matrix-method
#' kNNinterface,DataFrame-method kNNinterface,MultiAssayExperiment-method
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
#' \code{DataFrame} method or parameters that \code{\link[class]{knn}} can
#' accept.
#' @param classifierName Default: k Nearest Neighbours. Useful for automated
#' plot annotation by plotting functions within this package.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return A factor vector, the same as is returned by \code{\link[class]{knn}}.
#' @author Dario Strbenac
#' @examples
#' 
#'   if(require(class))
#'   {
#'     classes <- factor(rep(c("Healthy", "Disease"), each = 5), levels = c("Healthy", "Disease"))
#'     measurements <- matrix(c(rnorm(50, 10), rnorm(50, 5)), nrow = 10, byrow = TRUE)
#'     rownames(measurements) <- paste("Sample", 1:10)
#'     colnames(measurements) <- paste("mRNA", 1:10)
#'     
#'     # Train with 9 samples, test with one.
#'     kNNinterface(measurements[1:9, ], classes[1:9], measurements[10, , drop = FALSE])
#'   }
#'   
#' @usage NULL
#' @export
setGeneric("kNNinterface", function(measurementsTrain, ...) standardGeneric("kNNinterface"))

#' @rdname kNNinterface
#' @export
setMethod("kNNinterface", "matrix",
          function(measurementsTrain, classesTrain, measurementsTest, ...)
{
  kNNinterface(S4Vectors::DataFrame(measurementsTrain, check.names = FALSE),
               classesTrain,
               S4Vectors::DataFrame(measurementsTest, check.names = FALSE), ...)
})

#' @rdname kNNinterface
#' @export
setMethod("kNNinterface", "DataFrame", function(measurementsTrain, classesTrain, measurementsTest, ..., classifierName = "k Nearest Neighbours", verbose = 3)
{
  splitDataset <- .splitDataAndOutcome(measurementsTrain, classesTrain)
  classesTrain <- splitDataset[["outcome"]]
  trainingMatrix <- as.matrix(splitDataset[["measurements"]])
  measurementsTest <- measurementsTest[, colnames(measurementsTrain), drop = FALSE]
  
  if(!requireNamespace("class", quietly = TRUE))
    stop("The package 'class' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting k Nearest Neighbours classifier to data and predicting classes.")
  
  class::knn(as.matrix(measurementsTrain), as.matrix(measurementsTest), classesTrain, ...)
})

#' @rdname kNNinterface
#' @export
setMethod("kNNinterface", "MultiAssayExperiment",
function(measurementsTrain, measurementsTest, targets = names(measurementsTrain), classesTrain, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurementsTrain, targets, classesTrain)
  trainingTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["outcome"]]
  testingTable <- .MAEtoWideTable(measurementsTest, targets)
            
  .checkVariablesAndSame(trainingTable, testingTable)
  kNNinterface(trainingTable, classes, testingTable, ...)
})