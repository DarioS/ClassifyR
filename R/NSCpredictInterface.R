#' Interface for \code{pamr.predict} Function from \code{pamr} CRAN Package
#' 
#' Restructures variables from ClassifyR framework to be compatible with
#' \code{\link[pamr]{pamr.predict}} definition.
#' 
#' This function is an interface between the ClassifyR framework and
#' \code{\link[pamr]{pamr.predict}}.  It selects the highest threshold that
#' gives the minimum error rate in the training data.
#' 
#' @aliases NSCpredictInterface NSCpredictInterface,pamrtrained,matrix-method
#' NSCpredictInterface,pamrtrained,DataFrame-method
#' NSCpredictInterface,pamrtrained,MultiAssayExperiment-method
#' @param trained An object of class \code{pamrtrained}.
#' @param test An object of the same class as \code{measurements} with no
#' samples in common with \code{measurements} used in the training stage and
#' the same number of features as it.  Also, if a \code{DataFrame}, the
#' \code{class} column must be absent.
#' @param classes Either NULL or a character vector of length 1, specifying the
#' column name to remove.
#' @param targets If \code{test} is a \code{MultiAssayExperiment}, the names of
#' the data tables to be used. \code{"clinical"} is also a valid value and
#' specifies that numeric variables from the clinical data table will be used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method or optional settings that are passed to
#' \code{\link[pamr]{pamr.predict}}.
#' @param returnType Default: \code{"both"}. Either \code{"class"},
#' \code{"score"} or \code{"both"}.  Sets the return value from the prediction
#' to either a vector of class labels, score for a sample belonging to the
#' second class, as determined by the factor levels, or both labels and scores
#' in a \code{data.frame}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return Either a factor vector of predicted classes, a matrix of scores for
#' each class, or a table of both the class labels and class scores, depending
#' on the setting of \code{returnType}.
#' @author Dario Strbenac
#' @seealso \code{\link[pamr]{pamr.predict}} for the function that was
#' interfaced to.
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
#'     fit <- NSCtrainInterface(genesMatrix[, c(1:20, 26:45)], classes[c(1:20, 26:45)])
#'     NSCpredictInterface(fit, genesMatrix[, c(21:25, 46:50)])
#'   }
#' 
#' @export
setGeneric("NSCpredictInterface", function(trained, test, ...)
standardGeneric("NSCpredictInterface"))

setMethod("NSCpredictInterface", c("pamrtrained", "matrix"), function(trained, test, ...)
{
  NSCpredictInterface(trained, DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("NSCpredictInterface", c("pamrtrained", "DataFrame"), function(trained, test, classes = NULL, ..., returnType = c("both", "class", "score"), verbose = 3)
{
  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")
  returnType <- match.arg(returnType)
  
  if(!is.null(classes)) # Remove them.
  {
    splitDataset <- .splitDataAndClasses(test, classes) # Remove classes, if present.
    test <- splitDataset[["measurements"]]
  }
 
  minError <- min(trained[["errors"]])
  threshold <- trained[["threshold"]][max(which(trained[["errors"]] == minError))]

  test <- t(as.matrix(test))   
  classPredictions <- pamr::pamr.predict(trained, test, threshold, ...)
  classScores <- pamr::pamr.predict(trained, test, threshold, type = "posterior", ...)[, levels(trained[["y"]])]
  if(!is.matrix(classScores)) # Only one sample was predicted and pamr isn't consistent with return types.
    classScores <- t(classScores)
  
  if(verbose == 3)
    message("Nearest shrunken centroid predictions made.")
  
  switch(returnType, class = classPredictions, # Factor vector.
         score = classScores, # Numeric matrix.
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
})

setMethod("NSCpredictInterface", c("pamrtrained", "MultiAssayExperiment"), function(trained, test, targets = names(test), ...)
{
  test <- .MAEtoWideTable(test, targets)[["dataTable"]] # Remove any classes, if present.
  
  if(ncol(test) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    NSCpredictInterface(trained, test, ...)
})
