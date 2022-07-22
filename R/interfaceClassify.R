#' An Interface for PoiClaClu Package's Classify Function
#' 
#' More details of Poisson LDA are available in the documentation of
#' \code{\link[PoiClaClu]{Classify}}. Data tables which consist entirely of
#' non-integer data cannot be analysed.
#' 
#' @aliases classifyInterface classifyInterface,matrix-method
#' classifyInterface,DataFrame-method
#' classifyInterface,MultiAssayExperiment-method
#' @param countsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
#' If of type \code{\link{DataFrame}} or \code{\link{MultiAssayExperiment}}, the data set is subset
#' to only those features of type \code{numeric}.
#' @param classesTrain A vector of class labels of class \code{\link{factor}} of the
#' same length as the number of samples in \code{countsTrain} if it is a
#' \code{\link{matrix}} or a \code{\link{DataFrame}} or a character vector of length 1
#' containing the column name in \code{countsTrain} if it is a \code{\link{DataFrame}} or the
#' column name in \code{colData(countsTrain)} if \code{countsTrain} is a
#' \code{\link{MultiAssayExperiment}}. If a column name, that column will be
#' removed before training.
#' @param countsTest An object of the same class as \code{countsTrain} with no
#' samples in common with \code{countsTrain} and the same number of features
#' as it.
#' @param targets If \code{countsTrain} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"sampleInfo"} is also a valid value
#' and specifies that integer variables from the sample information table will be
#' used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method or parameters that \code{\link[PoiClaClu]{Classify}}
#' can accept.
#' @param returnType Default: \code{"both"}. Either \code{"class"},
#' \code{"score"} or \code{"both"}.  Sets the return value from the prediction
#' to either a vector of class labels, matrix of scores for each class, or both
#' labels and scores in a \code{data.frame}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return Either a factor vector of predicted classes, a matrix of scores for
#' each class, or a table of both the class labels and class scores, depending
#' on the setting of \code{returnType}.
#' @author Dario Strbenac
#' @examples
#' 
#'   if(require(PoiClaClu))
#'   {
#'     readCounts <- CountDataSet(n = 100, p = 1000, 2, 5, 0.1)
#'     # Rows are for features, columns are for samples.
#'     trainData <- readCounts[['x']]
#'     trainClasses <- factor(paste("Class", readCounts[['y']]))
#'     testData <- readCounts[['xte']]
#'     storage.mode(trainData) <- storage.mode(testData) <- "integer"
#'     classified <- classifyInterface(trainData, trainClasses, testData)
#'     
#'     setNames(table(paste("Class", readCounts[["yte"]]) == classified), c("Incorrect", "Correct"))
#'   }
#'   
#' @usage NULL
#' @export
setGeneric("classifyInterface", function(countsTrain, ...)
standardGeneric("classifyInterface"))

#' @rdname classifyInterface
#' @export
setMethod("classifyInterface", "matrix", # Matrix of integer measurements.
          function(countsTrain, classesTrain, countsTest, ...)
{
  classifyInterface(S4Vectors::DataFrame(countsTrain, check.names = FALSE),
                    classesTrain,
                    S4Vectors::DataFrame(countsTest, check.names = FALSE), ...)
})

# Sample information data or one of the other inputs, transformed.
#' @rdname classifyInterface
#' @export
setMethod("classifyInterface", "DataFrame", function(countsTrain, classesTrain, countsTest, ...,
                                returnType = c("both", "class", "score"), verbose = 3)
{
  if(!requireNamespace("PoiClaClu", quietly = TRUE))
    stop("The package 'PoiClaClu' could not be found. Please install it.")
  returnType <- match.arg(returnType)
  
  # Ensure that any non-integer variables are removed from the training and testing matrices.
  splitDataset <- .splitDataAndOutcome(countsTrain, classesTrain, restrict = "integer")
  classesTrain <- splitDataset[["outcome"]]
  trainingMatrix <- as.matrix(splitDataset[["measurements"]])
  isInteger <- sapply(countsTest, is.integer)
  testingMatrix <- as.matrix(countsTest[, isInteger, drop = FALSE])
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  
  if(verbose == 3)
    message("Fitting Poisson LDA classifier to training data and making predictions on test data.")

  predicted <- PoiClaClu::Classify(trainingMatrix, classesTrain, testingMatrix, ...)
  classPredictions <- predicted[["ytehat"]]
  classScores <- predicted[["discriminant"]]
  colnames(classScores) <- levels(classesTrain)
  switch(returnType, class = classPredictions, # Factor vector.
         score = classScores, # Numeric matrix.
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
})

#' @rdname classifyInterface
#' @export
setMethod("classifyInterface", "MultiAssayExperiment",
function(countsTrain, countsTest, targets = names(countsTrain), classesTrain, ...)
{
  tablesAndOutcome <- .MAEtoWideTable(countsTrain, targets, classesTrain, "integer")
  trainingMatrix <- tablesAndOutcome[["dataTable"]]
  classesTrain <- tablesAndOutcome[["outcome"]]
  testingMatrix <- .MAEtoWideTable(countsTest, targets, "integer")
            
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  classifyInterface(trainingMatrix, classesTrain, testingMatrix, ...)
})