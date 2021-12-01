setGeneric("calcExternalPerformance", function(actualClasses, predictedClasses, ...)
standardGeneric("calcExternalPerformance"))

setGeneric("calcCVperformance", function(result, ...)
standardGeneric("calcCVperformance"))

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
  levels(predictedClasses) <- levels(actualClasses)
  .calcPerformance(list(actualClasses), list(predictedClasses), performanceType = performanceType)[["values"]]
})

setMethod("calcCVperformance", "ClassifyResult",
          function(result, performanceType = c("Error", "Accuracy", "Balanced Error", "Balanced Accuracy",
                                               "Sample Error", "Sample Accuracy",
                                               "Micro Precision", "Micro Recall",
                                               "Micro F1", "Macro Precision",
                                               "Macro Recall", "Macro F1", "Matthews Correlation Coefficient"))
{
  performanceType <- match.arg(performanceType)
  if(length(levels(actualClasses)) > 2 && performanceType == "Matthews Correlation Coefficient")
    stop("Error: Matthews Correlation Coefficient specified but data set has more than 2 classes.")

  classLevels <- levels(actualClasses(result))
  samples <- factor(result@predictions[, "sample"], levels = sampleNames(result))
  predictedClasses <- factor(result@predictions[, "class"], levels = classLevels)
  actualClasses <- factor(actualClasses(result)[match(result@predictions[, "sample"], sampleNames(result))], levels = classLevels, ordered = TRUE)
  if(!performanceType %in% c("Sample Error", "Sample Accuracy"))
  {
    if("permutation" %in% colnames(result@predictions))
      grouping <- result@predictions[, "permutation"]
    else # A set of folds or all leave-k-out predictions.
      grouping <- rep(1, nrow(result@predictions))
  }
  performance <- .calcPerformance(actualClasses, predictedClasses, samples, performanceType, grouping)
  result@performance[[performance[["name"]]]] <- performance[["values"]]
  result
})

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
  } else { # Metrics for which true positives, true negatives, false positives, false negatives must be calculated.
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