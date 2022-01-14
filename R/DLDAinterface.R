#' An Interface for sparsediscrim Package's dlda Function
#' 
#' \code{DLDAtrainInterface} generates a trained diagonal LDA classifier and
#' \code{DLDApredictInterface} uses it to make predictions on a test data set.
#' 
#' If \code{measurements} is an object of class \code{MultiAssayExperiment},
#' the factor of sample classes must be stored in the DataFrame accessible by
#' the \code{colData} function with column name \code{"class"}.
#' 
#' @aliases DLDAtrainInterface DLDAtrainInterface,matrix-method
#' DLDAtrainInterface,DataFrame-method
#' DLDAtrainInterface,MultiAssayExperiment-method DLDApredictInterface
#' DLDApredictInterface,dlda,matrix-method
#' DLDApredictInterface,dlda,DataFrame-method
#' DLDApredictInterface,dlda,MultiAssayExperiment-method
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
#' @param model A fitted model as returned by \code{DLDAtrainInterface}.
#' @param test An object of the same class as \code{measurements} with no
#' samples in common with \code{measurements} and the same number of features
#' as it. Also, if a \code{DataFrame}, the \code{class} column must be absent.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"clinical"} is also a valid value
#' and specifies that integer variables from the clinical data table will be
#' used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method (e.g. \code{verbose}).
#' @param returnType Default: \code{"both"}. Either \code{"class"},
#' \code{"score"} or \code{"both"}.  Sets the return value from the prediction
#' to either a vector of class labels, matrix of scores for each class, or both
#' labels and scores in a \code{data.frame}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return For \code{DLDAtrainInterface}, a trained DLDA classifier.  For
#' \code{DLDApredictInterface}, either a factor vector of predicted classes, a
#' matrix of scores for each class, or a table of both the class labels and
#' class scores, depending on the setting of \code{returnType}.
#' @author Dario Strbenac
#' @examples
#' 
#'   # if(require(sparsediscrim)) Package currently removed from CRAN.
#'   #{
#'     # Genes 76 to 100 have differential expression.
#'     genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
#'     genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
#'                                       c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix))
#'     rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
#'     selected <- rownames(genesMatrix)[91:100]
#'     trainingSamples <- c(1:20, 26:45)
#'     testingSamples <- c(21:25, 46:50)
#'     
#'     classifier <- DLDAtrainInterface(genesMatrix[selected, trainingSamples],
#'                                      classes[trainingSamples])
#'     DLDApredictInterface(classifier, genesMatrix[selected, testingSamples])
#'   #}
#' 
#' @export
setGeneric("DLDAtrainInterface", function(measurements, ...)
standardGeneric("DLDAtrainInterface"))

setMethod("DLDAtrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  DLDAtrainInterface(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("DLDAtrainInterface", "DataFrame", function(measurements, classes, verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  trainingMatrix <- as.matrix(splitDataset[["measurements"]])
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  
  #if(!requireNamespace("sparsediscrim", quietly = TRUE))
    #stop("The package 'sparsediscrim' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting DLDA classifier to data.")
  
  # sparsediscrim::dlda(as.matrix(measurements), classes)
  .dlda(as.matrix(measurements), classes)
})

setMethod("DLDAtrainInterface", "MultiAssayExperiment",
function(measurements, targets = names(measurements), classes, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  
  if(ncol(measurements) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    DLDAtrainInterface(measurements, classes, ...)
})

setGeneric("DLDApredictInterface", function(model, test, ...)
standardGeneric("DLDApredictInterface"))

setMethod("DLDApredictInterface", c("dlda", "matrix"),
          function(model, test, ...)
{
  DLDApredictInterface(model, DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("DLDApredictInterface", c("dlda", "DataFrame"), function(model, test, returnType = c("both", "class", "score"), verbose = 3)
{
  isNumeric <- sapply(test, is.numeric)
  test <- test[, isNumeric, drop = FALSE]
  returnType <- match.arg(returnType)
  
  #if(!requireNamespace("sparsediscrim", quietly = TRUE)) # Removed from CRAN, sadly.
    #stop("The package 'sparsediscrim' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained DLDA classifier.")
  
  #predict(model, as.matrix(test))
  predictions <- .predict(model, as.matrix(test)) # Copy in utilities.R.
  
  switch(returnType, class = predictions[["class"]], # Factor vector.
                   score = predictions[["posterior"]][, model[["groups"]]], # Numeric matrix.
                   both = data.frame(class = predictions[["class"]], predictions[["posterior"]][, model[["groups"]]], check.names = FALSE))
})

setMethod("DLDApredictInterface", c("dlda", "MultiAssayExperiment"),
          function(model, test, targets = names(test), ...)
{
  tablesAndClasses <- .MAEtoWideTable(test, targets)
  test <- tablesAndClasses[["dataTable"]]
            
  if(ncol(test) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    DLDApredictInterface(model, test, ...)
})