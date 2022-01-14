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
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix}, the rows are features, and the columns are samples.  If of
#' type \code{\link{DataFrame}}, the data set is subset to only those features
#' of type \code{integer}.
#' @param classes Either a vector of class labels of class \code{\link{factor}}
#' of the same length as the number of samples in \code{measurements} or if the
#' measurements are of class \code{DataFrame} a character vector of length 1
#' containing the column name in \code{measurement} is also permitted. Not used
#' if \code{measurements} is a \code{MultiAssayExperiment} object.
#' @param test An object of the same class as \code{measurements} with no
#' samples in common with \code{measurements} and the same number of features
#' as it.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"clinical"} is also a valid value
#' and specifies that integer variables from the clinical data table will be
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
#' @return A factor vector, the same as is returned by
#' \code{\link[class]{knn}}.
#' @author Dario Strbenac
#' @examples
#' 
#'   if(require(class))
#'   {
#'     classes <- factor(rep(c("Healthy", "Disease"), each = 5),
#'                       levels = c("Healthy", "Disease"))
#'     measurements <- matrix(c(rnorm(50, 10), rnorm(50, 5)), ncol = 10)
#'     colnames(measurements) <- paste("Sample", 1:10)
#'     rownames(measurements) <- paste("mRNA", 1:10)
#'     
#'     kNNinterface(measurements[, 1:9], classes[1:9], measurements[, 10, drop = FALSE])
#'   }
#'   
#' @export
setGeneric("kNNinterface", function(measurements, ...) standardGeneric("kNNinterface"))

setMethod("kNNinterface", "matrix",
          function(measurements, classes, test, ...)
{
  kNNinterface(DataFrame(t(measurements), check.names = FALSE),
               classes,
               DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("kNNinterface", "DataFrame", function(measurements, classes, test, ..., classifierName = "k Nearest Neighbours", verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  trainingMatrix <- as.matrix(splitDataset[["measurements"]])
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  isNumeric <- sapply(test, is.numeric)
  test <- test[, isNumeric, drop = FALSE]  
  
  if(!requireNamespace("class", quietly = TRUE))
    stop("The package 'class' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting k Nearest Neighbours classifier to data and predicting classes.")
  
  class::knn(as.matrix(measurements), as.matrix(test), classes, ...)
})

setMethod("kNNinterface", "MultiAssayExperiment",
function(measurements, test, targets = names(measurements), classes, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes)
  trainingTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  testingTable <- .MAEtoWideTable(test, targets)
            
  .checkVariablesAndSame(trainingTable, testingTable)
  kNNinterface(trainingTable, classes, testingTable, ...)
})