#' Interface for \code{pamr.train} Function from \code{pamr} CRAN Package
#' 
#' Restructures variables from ClassifyR framework to be compatible with
#' \code{\link[pamr]{pamr.train}} definition.
#' 
#' This function is an interface between the ClassifyR framework and
#' \code{\link[pamr]{pamr.train}}.
#' 
#' @aliases NSCtrainInterface NSCtrainInterface,matrix-method
#' NSCtrainInterface,DataFrame-method
#' NSCtrainInterface,MultiAssayExperiment-method
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix}, the rows are features, and the columns are samples.
#' @param classes Either a vector of class labels of class \code{\link{factor}}
#' of the same length as the number of samples in \code{measurements} or if the
#' measurements are of class \code{DataFrame} a character vector of length 1
#' containing the column name in \code{measurement} is also permitted. Not used
#' if \code{measurements} is a \code{MultiAssayExperiment} object.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"clinical"} is also a valid value
#' and specifies that numeric variables from the clinical data table will be
#' used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method or extra arguments passed to
#' \code{\link[pamr]{pamr.train}}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return A list with elements as described in \code{\link[pamr]{pamr.train}}.
#' @author Dario Strbenac
#' @seealso \code{\link[pamr]{pamr.train}} for the function that was interfaced
#' to.
#' @examples
#' 
#'   if(require(pamr))
#'   {
#'     # Samples in one class with differential expression to other class.
#'     genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(100, 9, 1)))
#'     genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn)
#'                                  c(rnorm(75, 9, 1), rnorm(25, 14, 1))))
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     
#'     NSCtrainInterface(genesMatrix, classes)
#'   }
#' 
#' @export
setGeneric("NSCtrainInterface", function(measurements, ...)
standardGeneric("NSCtrainInterface"))

setMethod("NSCtrainInterface", "matrix", function(measurements, classes, ...)
{
  NSCtrainInterface(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("NSCtrainInterface", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, ..., verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")

  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")

  trainedModel <- pamr::pamr.train(list(x = t(as.matrix(measurements)), y = classes), ...)
  
  if(verbose == 3)
    message("Nearest shrunken centroid training completed.")
  
  trainedModel  
})

setMethod("NSCtrainInterface", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), classes, ...)
{ 
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  
  if(ncol(measurements) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    NSCtrainInterface(measurements, classes, ...)
})