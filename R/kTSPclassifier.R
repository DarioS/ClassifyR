#' Classification Using k Pairs of Features With Relative Differences Between
#' Classes
#' 
#' Each pair of features votes for a class based on whether the value of one
#' feature is less than the other feature.
#' 
#' If the voting is tied, the the class with the most samples in the training
#' set is voted for.
#' 
#' Because this method compares different features, they need to have
#' comparable measurements. For example, RNA-seq counts would be unsuitable
#' since these depend on the length of a feature, whereas F.P.K.M. values would
#' be suitable.
#' 
#' The \code{featurePairs} to use is recommended to be determined in
#' conjunction with \code{pairsDifferencesRanking}.
#' 
#' @aliases kTSPclassifier kTSPclassifier,matrix-method
#' kTSPclassifier,DataFrame-method kTSPclassifier,MultiAssayExperiment-method
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix}, the rows are features, and the columns are samples.
#' @param classes Either a vector of class labels of class \code{\link{factor}}
#' of the same length as the number of samples in \code{measurements} or if the
#' measurements are of class \code{DataFrame} a character vector of length 1
#' containing the column name in \code{measurement} is also permitted. Not used
#' if \code{measurements} is a \code{MultiAssayExperiment} object.
#' @param test An object of the same class as \code{measurements} with no
#' samples in common with \code{measurements} and the same number of features
#' as it.
#' @param featurePairs An object of class as \code{\link{Pairs}} containing the
#' pairs of features to determine whether the inequality of the first feature
#' being less than the second feature holds, indiciating evidence for the
#' second level of the \code{classes} factor.
#' @param difference Default: \code{"unweighted"}. Either \code{"unweighted"}
#' or \code{"weighted"}.  In weighted mode, the difference in densities is
#' summed over all features.  If unweighted mode, each feature's vote is worth
#' the same.
#' @param minDifference Default: 0. The minimum difference in densities for a
#' feature to be allowed to vote. Can be a vector of cutoffs. If no features
#' for a particular sample have a difference large enough, the class predicted
#' is simply the largest class.
#' @param target If \code{measurements} is a \code{MultiAssayExperiment}, the
#' name of the data table to be used. \code{"clinical"} is also a valid value
#' and specifies that integer variables from the clinical data table will be
#' used.
#' @param ... Unused variables by the methods for a \code{matrix} or a
#' \code{MultiAssayExperiment} passed to the \code{DataFrame} method which does
#' the classification.
#' @param returnType Default: \code{"both"}. Either \code{"class"},
#' \code{"score"} or \code{"both"}.  Sets the return value from the prediction
#' to either a vector of class labels, score for a sample belonging to the
#' second class, as determined by the factor levels, or both labels and scores
#' in a \code{data.frame}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return A vector or data frame of class prediction information, as long as
#' the number of samples in the test data.
#' @author Dario Strbenac
#' @seealso \code{\link{pairsDifferencesRanking}} for a function which could be
#' used to do feature ranking before the k-TSP classifier is run.
#' @examples
#' 
#'   # Difference in differences for features A and C between classes.                                           
#'   measurements <- matrix(c(9.9, 10.5, 10.1, 10.9, 11.0, 6.6, 7.7, 7.0, 8.1, 6.5,
#'                            8.5, 10.5, 12.5, 10.5, 9.5, 8.5, 10.5, 12.5, 10.5, 9.5,
#'                            6.6, 7.7, 7.0, 8.1, 6.5, 11.2, 11.0, 11.1, 11.4, 12.0,
#'                            8.1, 10.6, 7.4, 7.1, 10.4, 6.1, 7.3, 2.7, 11.0, 9.1,
#'                            round(rnorm(60, 8, 1), 1)), ncol = 10, byrow = TRUE)
#'   classes <- factor(rep(c("Good", "Poor"), each = 5))
#'                          
#'   rownames(measurements) <- LETTERS[1:10]
#'   colnames(measurements) <- names(classes) <- paste("Patient", 1:10)
#'   
#'   trainIndex <- c(1:4, 6:9)
#'   trainMatrix <- measurements[, trainIndex]
#'   testMatrix <- measurements[, c(5, 10)]
#'   
#'   featurePairs <- Pairs('A', 'C') # Could be ranked by pairsDifferencesRanking function and
#'                                   # selected internally by tuning of features selected.
#'   kTSPclassifier(trainMatrix, classes[trainIndex], testMatrix, featurePairs)
#' 
#' @export
setGeneric("kTSPclassifier", function(measurements, ...)
           standardGeneric("kTSPclassifier"))

setMethod("kTSPclassifier", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, test, featurePairs = NULL, ...)
{
  if(is.null(featurePairs))
    stop("No feature pairs provided but some must be.")
  if(!"Pairs" %in% class(featurePairs))
    stop("'featurePairs' must be of type Pairs.")            
            
  kTSPclassifier(DataFrame(t(measurements[, , drop = FALSE]), check.names = FALSE),
                 classes,
                 DataFrame(t(test[, , drop = FALSE]), check.names = FALSE), featurePairs, ...)
})

setMethod("kTSPclassifier", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, test, featurePairs = NULL,
                   difference = c("unweighted", "weighted"), minDifference = 0,
                   returnType = c("both", "class", "score"), verbose = 3)
{
  if(is.null(featurePairs))
    stop("No feature pairs provided but some must be.")
  if(!"Pairs" %in% class(featurePairs))
    stop("'featurePairs' must be of type Pairs.")            
          
  splitDataset <- .splitDataAndClasses(measurements, classes)
  trainingMatrix <- splitDataset[["measurements"]]
  isNumeric <- sapply(trainingMatrix, is.numeric)
  trainingMatrix <- as.matrix(trainingMatrix[, isNumeric, drop = FALSE])
  isNumeric <- sapply(test, is.numeric)
  testingMatrix <- as.matrix(test[, isNumeric, drop = FALSE])
  
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  
  difference <- match.arg(difference)
  returnType <- match.arg(returnType)
  
  classesSizes <- sapply(levels(classes), function(class) sum(classes == class))
  largerClass <- names(classesSizes)[which.max(classesSizes)[1]]
  secondClass <- classes == levels(classes)[2]
  
  if(verbose == 3)
    message("Determining inequalities of feature pairs.")

  # Order pairs so that first < second is the rule for predicting the second class, based on factor levels.
  # Effectively the classifier training.
  
  featurePairs <- do.call(c, lapply(featurePairs, function(pair)
  {
    isSmaller <- trainingMatrix[secondClass, S4Vectors::first(pair)] < trainingMatrix[secondClass, S4Vectors::second(pair)]
    if(sum(isSmaller) < length(isSmaller) / 2)
      Pairs(S4Vectors::second(pair), S4Vectors::first(pair))
    else
      pair
  }))

  testDataFrame <- data.frame(t(testingMatrix), check.names = FALSE)
  
  if(verbose == 3)
    message("Predicting sample classes using feature pair inequalities.")
  
  predictions <- do.call(rbind, lapply(testDataFrame, function(sampleValues)
  {
    names(sampleValues) <- rownames(testDataFrame)
    measureDifferences <- sampleValues[S4Vectors::second(featurePairs)] - sampleValues[S4Vectors::first(featurePairs)]
    useFeatures <- which(abs(measureDifferences) > minDifference)
    if(length(useFeatures) == 0) # No features have a large enough distance difference.
    {                            # Simply vote for the larger class.
      if(largerClass == levels(classes)[1])
      {
        class <- levels(classes)[1]
        score <- -1
      } else {
        class <- levels(classes)[2]
        score <- 1
      }
    } else { # One or more features are available to vote with.
      measureDifferences <- measureDifferences[useFeatures]
      if(difference == "unweighted")
      {
        # For being in second class.
        score <- sum(measureDifferences > 0)
            
        if(score > length(measureDifferences) / 2)
          class <- levels(classes)[2]
        else
          class <- levels(classes)[1]
            
      } else { # Each pair contributes a score for class prediction.
               # For being in second class.
        score <- sum(measureDifferences)

        # Sum of scores is tested for being positive or negative.
        class <- levels(classes)[(sum(measureDifferences) > 0) + 1]
      }
    }
    data.frame(class = factor(class, levels = levels(classes)), score = score, check.names = FALSE)
  }))

  switch(returnType, class = predictions[, "class"],
         score = predictions[, colnames(predictions) %in% levels(classes)],
         both = data.frame(class = predictions[, "class"], predictions[, colnames(predictions) %in% levels(classes)], check.names = FALSE)
  )
})

setMethod("kTSPclassifier", "MultiAssayExperiment", 
          function(measurements, test, target = names(measurements)[1], featurePairs = NULL, ...)
{
  if(is.null(featurePairs))
    stop("No feature pairs provided but some must be.")
  if(!"Pairs" %in% class(featurePairs))
    stop("'featurePairs' must be of type Pairs.")            
            
  tablesAndClasses <- .MAEtoWideTable(measurements, target)
  trainingMatrix <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  testingMatrix <- .MAEtoWideTable(test, target)
            
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  kTSPclassifier(trainingMatrix, classes, testingMatrix, featurePairs, ...)
})