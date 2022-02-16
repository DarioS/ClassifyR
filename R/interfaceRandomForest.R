################################################################################
#
# Train interface
#
################################################################################


#' An Interface for randomForest Package's randomForest Function
#' 
#' A random forest classifier builds multiple decision trees and uses the
#' predictions of the trees to determine a single prediction for each test
#' sample.
#' 
#' @aliases randomForestInterface randomForestTrainInterface
#' randomForestTrainInterface,matrix-method
#' randomForestTrainInterface,DataFrame-method
#' randomForestTrainInterface,MultiAssayExperiment-method
#' randomForestPredictInterface
#' randomForestPredictInterface,randomForest,matrix-method
#' randomForestPredictInterface,randomForest,DataFrame-method
#' randomForestPredictInterface,randomForest,MultiAssayExperiment-method
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
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
#'     genesMatrix <- sapply(1:100, function(sample) rnorm(25, 9, 0.3))
#'     genesMatrix <- rbind(genesMatrix, t(sapply(1:25, function(sample)
#'                                       c(rnorm(75, 9, 0.3), rnorm(25, 14, 0.3)))))
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix))
#'     rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
#'     trainingSamples <- c(1:20, 26:45)
#'     testingSamples <- c(21:25, 46:50)
#'     
#'     trained <- randomForestTrainInterface(genesMatrix[trainingSamples, ],
#'                                           classes[trainingSamples])
#'     predicted <- randomForestPredictInterface(trained, genesMatrix[testingSamples, ])
#'   }
#' 
#' @export
setGeneric("randomForestTrainInterface", function(measurementsTrain, ...)
standardGeneric("randomForestTrainInterface"))

#' @export
setGeneric("randomForestPredictInterface", function(models, measurementsTest, ...)
           standardGeneric("randomForestPredictInterface"))

setMethod("randomForestTrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurementsTrain, classesTrain, ...)
{
  randomForestTrainInterface(DataFrame(measurementsTrain, check.names = FALSE),
                             classesTrain, ...)
})

# Sample information data or one of the other inputs, transformed.
setMethod("randomForestTrainInterface", "DataFrame", function(measurementsTrain, classesTrain, ..., verbose = 3)
{
  splitDataset <- .splitDataAndOutcomes(measurementsTrain, classesTrain)

  if(!requireNamespace("randomForest", quietly = TRUE))
    stop("The package 'randomForest' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting random forest classifier to training data and making predictions on test
            data.")

  # Convert to base data.frame as randomForest doesn't understand DataFrame.
  randomForest::randomForest(as(splitDataset[["measurements"]], "data.frame"), splitDataset[["outcomes"]], keep.forest = TRUE, ...)
})

setMethod("randomForestTrainInterface", "MultiAssayExperiment",
function(measurementsTrain, targets = names(measurementsTrain), classesTrain, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurementsTrain, targets, classesTrain, restrict = NULL)
  measurementsTrain <- tablesAndClasses[["dataTable"]]
  classesTrain <- tablesAndClasses[["outcomes"]]
  
  randomForestTrainInterface(measurementsTrain, classesTrain, ...)
})



################################################################################
#
# Predict Interface
#
################################################################################



setGeneric("randomForestPredictInterface", function(forest, measurementsTest, ...)
           standardGeneric("randomForestPredictInterface"))

setMethod("randomForestPredictInterface", c("randomForest", "matrix"), function(forest, measurementsTest, ...)
{
  randomForestPredictInterface(forest, DataFrame(measurementsTest, check.names = FALSE), ...)
})

setMethod("randomForestPredictInterface", c("randomForest", "DataFrame"),
function(forest, measurementsTest, ..., returnType = c("both", "class", "score"), verbose = 3)
{
  returnType <- match.arg(returnType)
  if(verbose == 3)
    message("Predicting using random forest.")  
  
  classPredictions <- predict(forest, measurementsTest)
  classScores <- predict(forest, measurementsTest, type = "vote")[, forest[["classes"]], drop = FALSE]
  switch(returnType, class = classPredictions,
         score = classScores,
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
})

# One or more omics data sets, possibly with sample information data.
setMethod("randomForestPredictInterface", c("randomForest", "MultiAssayExperiment"),
          function(forest, measurementsTest, targets = names(measurementsTest), ...)
{
  testingTable <- .MAEtoWideTable(measurementsTest, targets)
  randomForestPredictInterface(models, testingTable, ...)
})




################################################################################
#
# Get selected features
#
################################################################################



#' Extract Vectors of Ranked and Selected Features From a Random Forest Object
#'
#' Provides a ranking of features based on the total decrease in node
#' impurities from splitting on the variable, averaged over all trees. Also
#' provides the selected features which are those that were used in at least
#' one tree of the forest.
#'
#'
#' @aliases forestFeatures forestFeatures,randomForest-method
#' @param forest A trained random forest which was created by
#' \code{\link{randomForest}}.
#' @return An \code{list} object. The first element is a vector or data frame
#' of features, ranked from best to worst using the Gini index. The second
#' element is a vector or data frame of features used in at least one tree.
#' @author Dario Strbenac
#' @examples
#'
#'       if(require(randomForest))
#'       {
#'         genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
#'         genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
#'                                           c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
#'         classes <- factor(rep(c("Poor", "Good"), each = 25))
#'         colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix))
#'         rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
#'         trainingSamples <- c(1:20, 26:45)
#'         testingSamples <- c(21:25, 46:50)
#'
#'         trained <- randomForestTrainInterface(genesMatrix[, trainingSamples],
#'                                               classes[trainingSamples], ntree = 10)
#'
#'         forestFeatures(trained)
#'       }
#'
#' @importFrom randomForest importance varUsed
#' @export
setGeneric("forestFeatures", function(forest, ...)
  standardGeneric("forestFeatures"))

setMethod("forestFeatures", "randomForest",
          function(forest)
          {
            inputFeatures <- rownames(randomForest::importance(forest))
            rankedFeatures <- inputFeatures[order(randomForest::importance(forest), decreasing = TRUE)]
            selectedFeatures <- inputFeatures[randomForest::varUsed(forest) > 0]
            selectedFeatures <- selectedFeatures[na.omit(match(rankedFeatures, selectedFeatures))]
            
            # Colon is a reserved symbol for separating data name and feature name, which is necessary for identifiability of MultiAssayExperiment features. It is not permitted in feature names.
            if(grepl(':', selectedFeatures[1]) == TRUE) # Convert to data.frame.
            {
              selectedFeatures <- do.call(rbind, strsplit(selectedFeatures, ':'))
              rankedFeatures <- do.call(rbind, strsplit(rankedFeatures, ':'))
              colnames(selectedFeatures) <- colnames(rankedFeatures) <- c("dataset", "feature")
            }
            list(rankedFeatures, selectedFeatures)
          })