#' Add Performance Calculations to a ClassifyResult Object or Calculate for a
#' Pair of Factor Vectors
#' 
#' If \code{calcExternalPerformance} is used, such as when having a vector of
#' known classes and a vector of predicted classes determined outside of the
#' ClassifyR package, a single metric value is calculated. If
#' \code{calcCVperformance} is used, annotates the results of calling
#' \code{\link{runTests}} with one of the user-specified performance measures.
#' 
#' All metrics except Matthews Correlation Coefficient are suitable for
#' evaluating classification scenarios with more than two classes and are
#' reimplementations of those available from
#' \href{https://software.intel.com/en-us/daal-programming-guide-details-40Intel
#' DAAL}.
#' 
#' If \code{\link{runTests}} was run in resampling mode, one performance
#' measure is produced for every resampling. If the leave-k-out mode was used,
#' then the predictions are concatenated, and one performance measure is
#' calculated for all classifications.
#' 
#' \code{"Balanced Error"} calculates the balanced error rate and is better
#' suited to class-imbalanced data sets than the ordinary error rate specified
#' by \code{"Error"}. \code{"Sample Error"} calculates the error rate of each
#' sample individually. This may help to identify which samples are
#' contributing the most to the overall error rate and check them for
#' confounding factors. Precision, recall and F1 score have micro and macro
#' summary versions. The macro versions are preferable because the metric will
#' not have a good score if there is substantial class imbalance and the
#' classifier predicts all samples as belonging to the majority class.
#' 
#' @aliases calcPerformance calcExternalPerformance calcCVperformance
#' calcExternalPerformance,factor,factor-method
#' calcCVperformance,ClassifyResult-method
#' @param result An object of class \code{\link{ClassifyResult}}.
#' @param performanceType A character vector of length 1. Default:
#' \code{"Balanced Error"}.\cr Must be one of the following options:\cr
#' \itemize{ \item\code{"Error"}: Ordinary error rate.  \item\code{"Accuracy"}:
#' Ordinary accuracy.  \item\code{"Balanced Error"}: Balanced error rate.
#' \item\code{"Balanced Accuracy"}: Balanced accuracy.  \item\code{"Sample
#' Error"}: Error rate for each sample in the data set.  \item\code{"Sample
#' Accuracy"}: Accuracy for each sample in the data set.
#' 
#' \item\code{"Micro Precision"}: Sum of the number of correct predictions in
#' each class, divided by the sum of number of samples in each class.
#' \item\code{"Micro Recall"}: Sum of the number of correct predictions in each
#' class, divided by the sum of number of samples predicted as belonging to
#' each class.  \item\code{"Micro F1"}: F1 score obtained by calculating the
#' harmonic mean of micro precision and micro recall.  \item\code{"Macro
#' Precision"}: Sum of the ratios of the number of correct predictions in each
#' class to the number of samples in each class, divided by the number of
#' classes.  \item\code{"Macro Recall"}: Sum of the ratios of the number of
#' correct predictions in each class to the number of samples predicted to be
#' in each class, divided by the number of classes.  \item\code{"Macro F1"}: F1
#' score obtained by calculating the harmonic mean of macro precision and macro
#' recall.  \item\code{"Matthews Correlation Coefficient"}: Matthews
#' Correlation Coefficient (MCC). A score between -1 and 1 indicating how
#' concordant the predicted classes are to the actual classes. Only defined if
#' there are two classes.  }
#' 
#' @param actualClasses A factor vector specifying each sample's correct class.
#' @param predictedClasses A factor vector of the same length as
#' \code{actualClasses} specifying each sample's predicted class.
#' 
#' @return If \code{calcCVperformance} was run, an updated
#' \code{\linkS4class{ClassifyResult}} object, with new metric values in the
#' \code{performance} slot. If \code{calcExternalPerformance} was run, the
#' performance metric value itself.
#' 
#' @author Dario Strbenac
#' @examples
#' 
#'   predictTable <- data.frame(sample = paste("A", 1:10, sep = ''),
#'                              class = factor(sample(LETTERS[1:2], 50, replace = TRUE)))
#'   actual <- factor(sample(LETTERS[1:2], 10, replace = TRUE))                             
#'   result <- ClassifyResult(DataFrame(),
#'                            paste("A", 1:10, sep = ''), paste("Gene", 1:50, sep = ''),
#'                            list(1:50, 1:50), list(1:5, 6:15), list(function(oracle){}), NULL,
#'                            predictTable, actual)
#'   result <- calcCVperformance(result) 
#'   performance(result)
#' 
#' @include classes.R
#' @export
setGeneric("calcExternalPerformance", function(actualClasses, predictedClasses, ...)
standardGeneric("calcExternalPerformance"))

setMethod("calcExternalPerformance", c("factor", "factor"),
          function(actualClasses, predictedClasses,
                   performanceType = c("Error", "Accuracy", "Balanced Error", "Balanced Accuracy",
                                       "Sample Error", "Sample Accuracy",
                                       "Micro Precision", "Micro Recall",
                                       "Micro F1", "Macro Precision",
                                       "Macro Recall", "Macro F1", "Matthews Correlation Coefficient"))
{
  performanceType <- match.arg(performanceType)
  if(length(levels(actualClasses)) > 2 && performanceType == "Matthews Correlation Coefficient")
    stop("Error: Matthews Correlation Coefficient specified but data set has more than 2 classes.")
  if(is(predictedClasses, "factor")) levels(predictedClasses) <- levels(actualClasses)
  .calcPerformance(list(actualClasses), list(predictedClasses), performanceType = performanceType)[["values"]]
})

setMethod("calcExternalPerformance", c("Surv", "numeric"),
          function(actualClasses, predictedClasses,
                   performanceType = c("C index"))
          {
            performanceType <- match.arg(performanceType)
            .calcPerformance(actualClasses, predictedClasses, performanceType = performanceType)[["values"]]
          })


#' @export
setGeneric("calcCVperformance", function(result, ...)
    standardGeneric("calcCVperformance"))


setMethod("calcCVperformance", "ClassifyResult",
          function(result, performanceType = c("Error", "Accuracy", "Balanced Error", "Balanced Accuracy",
                                               "Sample Error", "Sample Accuracy",
                                               "Micro Precision", "Micro Recall",
                                               "Micro F1", "Macro Precision",
                                               "Macro Recall", "Macro F1", "Matthews Correlation Coefficient", 
                                               "C index"))
{
  performanceType <- match.arg(performanceType)
  
  
  ### Group by permutation
  if(!performanceType %in% c("Sample Error", "Sample Accuracy"))
  {
    if("permutation" %in% colnames(result@predictions))
      grouping <- result@predictions[, "permutation"]
    else # A set of folds or all leave-k-out predictions.
      grouping <- rep(1, nrow(result@predictions))
  }
  
  ### Performance for survival data
  if(performanceType == "C index"){
    samples <- factor(result@predictions[, "sample"], levels = sampleNames(result))
    performance <- .calcPerformance(actualClasses = actualClasses(result)[match(result@predictions[, "sample"], sampleNames(result))],
                                    predictedClasses = result@predictions[, "class"], 
                                    samples = samples,
                                    performanceType = performanceType, 
                                    grouping = grouping)
    result@performance[[performance[["name"]]]] <- performance[["values"]]
    return(result)
  }
  
  
  ### Performance for data with classes
  if(length(levels(actualClasses)) > 2 && performanceType == "Matthews Correlation Coefficient")
    stop("Error: Matthews Correlation Coefficient specified but data set has more than 2 classes.")

  classLevels <- levels(actualClasses(result))
  samples <- factor(result@predictions[, "sample"], levels = sampleNames(result))
  predictedClasses <- factor(result@predictions[, "class"], levels = classLevels)
  actualClasses <- factor(actualClasses(result)[match(result@predictions[, "sample"], sampleNames(result))], levels = classLevels, ordered = TRUE)
  performance <- .calcPerformance(actualClasses, predictedClasses, samples, performanceType, grouping)
  result@performance[[performance[["name"]]]] <- performance[["values"]]
  result
})

#' @importFrom survival concordance
.calcPerformance <- function(actualClasses, predictedClasses, samples = NA, performanceType, grouping = NULL)
{
  if(performanceType %in% c("Sample Error", "Sample Accuracy"))
  {
    resultTable <- data.frame(sample = samples,
                              actual = actualClasses,
                              predicted = predictedClasses)
    allIDs <- levels(resultTable[, "sample"])
    sampleMetricValues <- sapply(allIDs, function(sampleID)
    {
      sampleResult <- subset(resultTable, sample == sampleID)
      if(nrow(sampleResult) == 0)
        return(NA)
      if(performanceType == "Sample Error")
        sum(as.character(sampleResult[, "predicted"]) != as.character(sampleResult[, "actual"]))
      else
        sum(as.character(sampleResult[, "predicted"]) == as.character(sampleResult[, "actual"]))
    })
    performanceValues <- as.numeric(sampleMetricValues / table(factor(resultTable[, "sample"], levels = allIDs)))
    names(performanceValues) <- allIDs
    return(list(name = performanceType, values = performanceValues))
  }
  

  
  if(!is.null(grouping))
  {
    actualClasses <- split(actualClasses, grouping)
    predictedClasses <- split(predictedClasses, grouping)
  }
  
  if(!is(actualClasses,"list")) actualClasses <- list(actualClasses)
  if(!is(predictedClasses,"list")) predictedClasses <- list(predictedClasses)
  

  if(performanceType %in% c("Accuracy", "Error")) {
    performanceValues <- unlist(mapply(function(iterationClasses, iterationPredictions)
    {
      # Columns are predicted classes, rows are actual classes.
      confusionMatrix <- table(iterationClasses, iterationPredictions)
      totalPredictions <- sum(confusionMatrix)
      correctPredictions <- sum(diag(confusionMatrix))
      diag(confusionMatrix) <- 0
      wrongPredictions <- sum(confusionMatrix)
      if(performanceType == "Accuracy")
        correctPredictions / totalPredictions
      else # It is "error".
        wrongPredictions / totalPredictions
    }, actualClasses, predictedClasses, SIMPLIFY = FALSE))
  } else if(performanceType %in% c("Balanced Accuracy", "Balanced Error")) {
    performanceValues <- unlist(mapply(function(iterationClasses, iterationPredictions)
    {
      # Columns are predicted classes, rows are actual classes.
      confusionMatrix <- table(iterationClasses, iterationPredictions)
      classSizes <- rowSums(confusionMatrix)
      classErrors <- classSizes - diag(confusionMatrix)
      if(performanceType == "Balanced Accuracy")
        mean(diag(confusionMatrix) / classSizes)
      else
        mean(classErrors / classSizes)
    }, actualClasses, predictedClasses, SIMPLIFY = FALSE))
  } else if(performanceType %in% c("C index")){
    performanceValues <- unlist(mapply(function(x,y){
      y <- -y
      survival::concordance(x ~ y)$concordance
    },actualClasses, predictedClasses, SIMPLIFY = FALSE))

    }else { # Metrics for which true positives, true negatives, false positives, false negatives must be calculated.
    performanceValues <- unlist(mapply(function(iterationClasses, iterationPredictions)
    {
      confusionMatrix <- table(iterationClasses, iterationPredictions)
      truePositives <- diag(confusionMatrix)
      falsePositives <- colSums(confusionMatrix) - truePositives
      falseNegatives <- rowSums(confusionMatrix) - truePositives
      trueNegatives <- sum(truePositives) - truePositives
      
      if(performanceType == "Average Accuracy")
      {
        return(sum((truePositives + trueNegatives) / (truePositives + trueNegatives + falsePositives + falseNegatives)) / nrow(confusionMatrix))
      }
      if(performanceType %in% c("Micro Precision", "Micro F1"))
      {
        microP <- sum(truePositives) / sum(truePositives + falsePositives)
        if(performanceType == "Micro Precision") return(microP)
      }
      if(performanceType %in% c("Micro Recall", "Micro F1"))
      {
        microR <- sum(truePositives) / sum(truePositives + falseNegatives)
        if(performanceType == "Micro Recall") return(microR)
      }
      if(performanceType == "Micro F1")
      {
        return(2 * microP * microR / (microP + microR))
      }
      if(performanceType %in% c("Macro Precision", "Macro F1"))
      {
        macroP <- sum(truePositives / (truePositives + falsePositives)) / nrow(confusionMatrix)
        if(performanceType == "Macro Precision") return(macroP)
      }
      if(performanceType %in% c("Macro Recall", "Macro F1"))
      {
        macroR <- sum(truePositives / (truePositives + falseNegatives)) / nrow(confusionMatrix)
        if(performanceType == "Macro Recall") return(macroR)
      }
      if(performanceType == "Macro F1")
      {
        return(2 * macroP * macroR / (macroP + macroR))
      }
      if(performanceType == "Matthews Correlation Coefficient")
      {
        return(unname((truePositives[2] * trueNegatives[2] - falsePositives[2] * falseNegatives[2]) / sqrt((truePositives[2] + falsePositives[2]) * (truePositives[2] + falseNegatives[2]) * (trueNegatives[2] + falsePositives[2]) * (trueNegatives[2] + falseNegatives[2]))))
      }
      
    }, actualClasses, predictedClasses, SIMPLIFY = FALSE))
  }

  list(name = performanceType, values = performanceValues)
}