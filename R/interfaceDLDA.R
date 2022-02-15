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
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
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
#' @param model A fitted model as returned by \code{DLDAtrainInterface}.
#' @param measurementsTest An object of the same class as \code{measurementsTrain} with no
#' samples in common with \code{measurementsTrain} and the same number of features
#' as it. Also, if a \code{DataFrame}, the \code{class} column must be absent.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"sampleInfo"} is also a valid value
#' and specifies that integer variables from the sample information data table will be
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
#' @include classes.R
#' @export
setGeneric("DLDAtrainInterface", function(measurementsTrain, ...) standardGeneric("DLDAtrainInterface"))

setMethod("DLDAtrainInterface", "matrix", function(measurementsTrain, classesTrain, ...) # Matrix of numeric measurements.
{
  DLDAtrainInterface(DataFrame(measurementsTrain, check.names = FALSE), classesTrain, ...)
})

setMethod("DLDAtrainInterface", "DataFrame", function(measurementsTrain, classesTrain, verbose = 3)
{
  splitDataset <- .splitDataAndOutcomes(measurementsTrain, classesTrain)
  trainingMatrix <- as.matrix(splitDataset[["measurements"]]) # DLDA demands matrix input type.
  classesTrain <- splitDataset[["outcomes"]]
  
  #if(!requireNamespace("sparsediscrim", quietly = TRUE))
  #stop("The package 'sparsediscrim' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting DLDA classifier to data.")
  
  # sparsediscrim::dlda(as.matrix(measurements), classes)
  .dlda(as.matrix(measurementsTrain), classesTrain)
})

setMethod("DLDAtrainInterface", "MultiAssayExperiment", function(measurementsTrain, targets = names(measurementsTrain), classesTrain, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurementsTrain, targets, classesTrain)
  measurementsTrain <- tablesAndClasses[["dataTable"]]
  classesTrain <- tablesAndClasses[["outcomes"]]
  
  if(ncol(measurementsTrain) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    DLDAtrainInterface(measurementsTrain, classesTrain, ...)
})


#' @export
setGeneric("DLDApredictInterface", function(model, measurementsTest, ...) standardGeneric("DLDApredictInterface"))

setMethod("DLDApredictInterface", c("dlda", "matrix"), function(model, measurementsTest, ...)
{
  DLDApredictInterface(model, DataFrame(measurementsTest, check.names = FALSE), ...)
})

setMethod("DLDApredictInterface", c("dlda", "DataFrame"), function(model, measurementsTest, returnType = c("both", "class", "score"), verbose = 3)
{
  isNumeric <- sapply(measurementsTest, is.numeric)
  measurementsTest <- measurementsTest[, isNumeric, drop = FALSE]
  returnType <- match.arg(returnType)
  
  #if(!requireNamespace("sparsediscrim", quietly = TRUE)) # Removed from CRAN, sadly.
  #stop("The package 'sparsediscrim' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained DLDA classifier.")
  
  #predict(model, as.matrix(test))
  predictions <- .predict(model, as.matrix(measurementsTest)) # Copy in utilities.R.
  
  switch(returnType, class = predictions[["class"]], # Factor vector.
         score = predictions[["posterior"]][, model[["groups"]]], # Numeric matrix.
         both = data.frame(class = predictions[["class"]], predictions[["posterior"]][, model[["groups"]]], check.names = FALSE))
})

setMethod("DLDApredictInterface", c("dlda", "MultiAssayExperiment"), function(model, measurementsTest, targets = names(measurementsTest), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurementsTest, targets)
  measurementsTest <- tablesAndClasses[["dataTable"]]
            
  if(ncol(measurementsTest) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    DLDApredictInterface(model, measurementsTest, ...)
})