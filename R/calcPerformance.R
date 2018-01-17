setGeneric("calcExternalPerformance", function(actualClasses, predictedClasses, ...)
{standardGeneric("calcExternalPerformance")})

setGeneric("calcCVperformance", function(result, ...)
{standardGeneric("calcCVperformance")})

setMethod("calcExternalPerformance", c("factor", "factor"),
          function(actualClasses, predictedClasses,
                   performanceType = c("error", "balanced error", "average accuracy",
                                       "micro precision", "micro recall",
                                       "micro F1", "macro precision",
                                       "macro recall", "macro F1"))
{
  performanceType <- match.arg(performanceType)
  .calcPerformance(list(actualClasses), list(predictedClasses), performanceType = performanceType)[["values"]]
})

setMethod("calcCVperformance", c("ClassifyResult"),
          function(result, performanceType = c("error", "balanced error", "sample error",
                                               "sample accuracy", "average accuracy",
                                               "micro precision", "micro recall",
                                               "micro F1", "macro precision",
                                               "macro recall", "macro F1"))
{
  performanceType <- match.arg(performanceType)
  classLevels <- levels(actualClasses(result))
  samples <- lapply(result@predictions, function(sample) factor(sample[, "sample"], levels = sampleNames(result)))
  predictedClasses <- lapply(result@predictions, function(sample) sample[, "label"])
  actualClasses <- lapply(result@predictions, function(sample)
                   factor(actualClasses(result)[sample[, "sample"]], levels = classLevels, ordered = TRUE))

  performance <- .calcPerformance(actualClasses, predictedClasses, samples, performanceType)
  
  result@performance[[performance[["name"]]]] <- performance[["values"]]
  result
})

.calcPerformance <- function(actualClasses, predictedClasses, samples = NA, performanceType)
{
  if(performanceType %in% c("sample error", "sample accuracy"))
  {
    resultTable <- data.frame(sample = unlist(samples),
                              actual = unlist(actualClasses),
                              predicted = unlist(predictedClasses))
    allIDs <- levels(resultTable[, "sample"])
    sampleMetricValues <- sapply(allIDs, function(sampleID)
    {
      sampleResult <- subset(resultTable, sample == sampleID)
      if(nrow(sampleResult) == 0)
        return(NA)
      if(performanceType == "sample error")
        sum(sampleResult[, "label"] != sampleResult[, "actual"])
      else
        sum(sampleResult[, "label"] == sampleResult[, "actual"])
    })
    performanceValues <- as.numeric(sampleMetricValues / table(factor(resultTable[, "sample"], levels = allIDs)))
    names(performanceValues) <- sampleNames(result)
    performanceName <- ifelse(performanceType == "sample error", "Sample-wise Error Rate", "Sample-wise Accuracy")
  } else if(performanceType == "error") {
    performanceValues <- unlist(mapply(function(iterationClasses, iterationPredictions)
    {
      # Columns are predicted classes, rows are actual classes.
      confusionMatrix <- table(iterationClasses, iterationPredictions)
      totalPredictions <- sum(confusionMatrix)
      diag(confusionMatrix) <- 0
      wrongPredictions <- sum(confusionMatrix)
      wrongPredictions / totalPredictions
    }, actualClasses, predictedClasses, SIMPLIFY = FALSE))
    performanceName <- "Error Rate"    
  } else if(performanceType == "balanced error") {
    performanceValues <- unlist(mapply(function(iterationClasses, iterationPredictions)
    {
      # Columns are predicted classes, rows are actual classes.
      confusionMatrix <- table(iterationClasses, iterationPredictions)
      classSizes <- rowSums(confusionMatrix)
      classErrors <- classSizes - diag(confusionMatrix)
      mean(classErrors / classSizes)
    }, actualClasses, predictedClasses, SIMPLIFY = FALSE))
    performanceName <- "Balanced Error Rate"
  } else { # Metrics for which true positives, true negatives, false positives, false negatives must be calculated.
    performanceName <- switch(performanceType, `average accuracy` = "Average Accuracy",
                              `micro precision` = "Micro Precision",
                              `micro recall` = "Micro Recall",
                              `micro F1` = "Micro F1 Score",
                              `macro precision` = "Macro Precision",
                              `macro recall` = "Macro Recall",
                              `macro F1` = "Macro F1 Score")
    performanceValues <- unlist(mapply(function(iterationClasses, iterationPredictions)
    {
      confusionMatrix <- table(iterationClasses, iterationPredictions)
      truePositives <- diag(confusionMatrix)
      falsePositives <- colSums(confusionMatrix) - truePositives
      falseNegatives <- rowSums(confusionMatrix) - truePositives
      trueNegatives <- sum(truePositives) - truePositives
      
      if(performanceType == "average accuracy")
      {
        return(sum((truePositives + trueNegatives) / (truePositives + trueNegatives + falsePositives + falseNegatives)) / nrow(confusionMatrix))
      }
      if(performanceType %in% c("micro precision", "micro F1"))
      {
        microP <- sum(truePositives) / sum(truePositives + falsePositives)
        if(performanceType == "micro precision") return(microP)
      }
      if(performanceType %in% c("micro recall", "micro F1"))
      {
        microR <- sum(truePositives) / sum(truePositives + falseNegatives)
        if(performanceType == "micro recall") return(microR)
      }
      if(performanceType == "micro F1")
      {
        return(2 * microP * microR / (microP + microR))
      }
      if(performanceType %in% c("macro precision", "macro F1"))
      {
        macroP <- sum(truePositives / (truePositives + falsePositives)) / nrow(confusionMatrix)
        if(performanceType == "macro precision") return(macroP)
      }
      if(performanceType %in% c("macro recall", "macro F1"))
      {
        marcoR <- sum(truePositives / (truePositives + falseNegatives)) / nrow(confusionMatrix)
        if(performanceType == "macro recall") return(macroR)
      }
      if(performanceType == "macro F1")
      {
        return(2 * macroP * marcoR / (macroP + marcoR))
      }
    }, actualClasses, predictedClasses, SIMPLIFY = FALSE))
  }
  list(name = performanceName, values = performanceValues)
}