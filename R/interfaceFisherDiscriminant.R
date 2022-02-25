#' Classification Using Fisher's LDA
#' 
#' Finds the decision boundary using the training set, and gives predictions
#' for the test set. Unlike ordinary LDA, Fisher's version does not have assumptions
#' about the normality of the features. Data tables which consist entirely of non-numeric
#' data cannot be analysed.
#' 
#' @aliases fisherDiscriminant fisherDiscriminant,matrix-method
#' fisherDiscriminant,DataFrame-method
#' fisherDiscriminant,MultiAssayExperiment-method
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
#' and specifies that numeric variables from the sample information data table will be
#' used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method.
#' @param returnType Default: \code{"both"}. Either \code{"class"},
#' \code{"score"}, or \code{"both"}.  Sets the return value from the prediction
#' to either a vector of class labels, score for a sample belonging to the
#' second class, as determined by the factor levels, or both labels and scores
#' in a \code{data.frame}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give. This function only prints progress messages if
#' the value is 3.
#' @return A vector or \code{data.frame} of class prediction information, as
#' long as the number of samples in the test data.
#' @author Dario Strbenac
#' @examples
#' 
#'   trainMatrix <- matrix(rnorm(1000, 8, 2), nrow = 10)
#'   classes <- factor(rep(c("Poor", "Good"), each = 5))
#'   
#'   # Make first 30 genes increased in value for poor samples.
#'   trainMatrix[1:5, 1:30] <- trainMatrix[1:5, 1:30] + 5
#'   
#'   testMatrix <- matrix(rnorm(1000, 8, 2), nrow = 10)
#'   
#'   # Make first 30 genes increased in value for sixth to tenth samples.
#'   testMatrix[6:10, 1:30] <- testMatrix[6:10, 1:30] + 5
#'   
#'   fisherDiscriminant(trainMatrix, classes, testMatrix)
#' 
#' @export
setGeneric("fisherDiscriminant", function(measurementsTrain, ...) standardGeneric("fisherDiscriminant"))

setMethod("fisherDiscriminant", "matrix", function(measurementsTrain, classesTrain, measurementsTest, ...) # Matrix of numeric measurements.
{
  fisherDiscriminant(DataFrame(measurementsTrain[, , drop = FALSE], check.names = FALSE),
                     classesTrain,
                     DataFrame(measurementsTest[, , drop = FALSE], check.names = FALSE), ...)
})

setMethod("fisherDiscriminant", "DataFrame", # Sample information data or one of the other inputs, transformed.
          function(measurementsTrain, classesTrain, measurementsTest, returnType = c("both", "class", "score"), verbose = 3)
{
  splitDataset <- .splitDataAndOutcomes(measurementsTrain, classesTrain)
  trainingMatrix <- as.matrix(splitDataset[["measurements"]])
  classesTrain <- splitDataset[["outcomes"]]
  isNumeric <- sapply(measurementsTest, is.numeric)
  testingMatrix <- as.matrix(measurementsTest[, isNumeric, drop = FALSE])
            
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  returnType <- match.arg(returnType)

  oneClassTraining <- which(classesTrain == levels(classesTrain)[1])
  otherClassTraining <- which(classesTrain == levels(classesTrain)[2])
  varOneClass <- apply(measurements[oneClassTraining, ], 2, var)
  varOtherClass <- apply(measurements[otherClassTraining, ], 2, var)
  varAll <- ((length(varOneClass) - 1) * varOneClass + (length(varOtherClass) - 1)
             * varOtherClass) / (length(oneClassTraining) + length(otherClassTraining) - 2)
  aT <- (apply(measurements[oneClassTraining, ], 2, mean) - apply(measurements[otherClassTraining, ], 2, mean)) / varAll
  criticalValue <- 0.5 * aT %*% as.matrix(apply(measurements[oneClassTraining, ], 2, mean) +
                                          apply(measurements[otherClassTraining, ], 2, mean))
  
  if(verbose == 3)
    message("Critical value calculated.")
  
  classesPredicted <- factor(apply(measurementsTest, 1, function(testSample)
  {
    if(aT %*% as.matrix(testSample) >= criticalValue)
      levels(classesTrain)[1]
    else
      levels(classesTrain)[2]
  }), levels = levels(classesTrain))
  scores <- apply(measurementsTest, 1, function(testSample) -1 * (aT %*% as.matrix(testSample))) # In reference to the second level of 'classes'. 
  
  switch(returnType, class = classesPredicted,
                     score = scores,
                     both = data.frame(class = classesPredicted, score = scores, check.names = FALSE))  
})

# One or more omics data sets, possibly with sample information data.
setMethod("fisherDiscriminant", "MultiAssayExperiment", function(measurementsTrain, measurementsTest, targets = names(measurementsTrain), classesTrain, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classesTrain)
  trainingMatrix <- tablesAndClasses[["dataTable"]]
  classesTrain <- tablesAndClasses[["outcomes"]]
  testingMatrix <- .MAEtoWideTable(measurementsTest, targets)
  
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  fisherDiscriminant(trainingMatrix, classesTrain, testingMatrix, ...)
})