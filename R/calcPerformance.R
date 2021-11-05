setGeneric("calcExternalPerformance", function(actualClasses, predictedClasses, ...)
standardGeneric("calcExternalPerformance"))

setGeneric("calcCVperformance", function(result, ...)
standardGeneric("calcCVperformance"))

setMethod("calcExternalPerformance", c("factor", "factor"),
          function(actualClasses, predictedClasses,
                   performanceType = c("error", "accuracy", "balanced error", "balanced accuracy",
                                       "sample error", "sample accuracy",
                                       "micro precision", "micro recall",
                                       "micro F1", "macro precision",
                                       "macro recall", "macro F1", "matthews"))
{
  performanceType <- match.arg(performanceType)
  if(length(levels(actualClasses)) > 2 && performanceType == "matthews")
    stop("Error: Matthews Correlation Coefficient specified but data set has more than 2 classes.")
  levels(predictedClasses) <- levels(actualClasses)
  .calcPerformance(list(actualClasses), list(predictedClasses), performanceType = performanceType)[["values"]]
})

setMethod("calcCVperformance", c("ClassifyResult"),
          function(result, performanceType = c("error", "accuracy", "balanced error", "balanced accuracy",
                                               "sample error", "sample accuracy",
                                               "micro precision", "micro recall",
                                               "micro F1", "macro precision",
                                               "macro recall", "macro F1", "matthews"))
{
  performanceType <- match.arg(performanceType)
  if(length(levels(actualClasses)) > 2 && performanceType == "matthews")
    stop("Error: Matthews Correlation Coefficient specified but data set has more than 2 classes.")
  
  classLevels <- levels(actualClasses(result))
  samples <- lapply(result@predictions, function(sample) factor(sample[, "sample"], levels = sampleNames(result)))
  predictedClasses <- lapply(result@predictions, function(sample) factor(sample[, "class"], levels = classLevels))
  actualClasses <- lapply(result@predictions, function(sample)
                   factor(actualClasses(result)[match(sample[, "sample"], sampleNames(result))], levels = classLevels, ordered = TRUE))

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
        sum(sampleResult[, "predicted"] != sampleResult[, "actual"])
      else
        sum(sampleResult[, "predicted"] == sampleResult[, "actual"])
    })
    performanceValues <- as.numeric(sampleMetricValues / table(factor(resultTable[, "sample"], levels = allIDs)))
    names(performanceValues) <- allIDs
    performanceName <- ifelse(performanceType == "sample error", "Sample-wise Error Rate", "Sample-wise Accuracy")
  } else if(performanceType %in% c("accuracy", "error")) {
    performanceValues <- unlist(mapply(function(iterationClasses, iterationPredictions)
    {
      # Columns are predicted classes, rows are actual classes.
      confusionMatrix <- table(iterationClasses, iterationPredictions)
      totalPredictions <- sum(confusionMatrix)
      correctPredictions <- sum(diag(confusionMatrix))
      diag(confusionMatrix) <- 0
      wrongPredictions <- sum(confusionMatrix)
      if(performanceType == "accuracy")
        correctPredictions / totalPredictions
      else # It is "error".
        wrongPredictions / totalPredictions
    }, actualClasses, predictedClasses, SIMPLIFY = FALSE))
    if(performanceType == "accuracy")
      performanceName <- "Accuracy"
    else # It is "error".
      performanceName <- "Error Rate"
  } else if(performanceType %in% c("balanced accuracy", "balanced error")) {
    performanceValues <- unlist(mapply(function(iterationClasses, iterationPredictions)
    {
      # Columns are predicted classes, rows are actual classes.
      confusionMatrix <- table(iterationClasses, iterationPredictions)
      classSizes <- rowSums(confusionMatrix)
      classErrors <- classSizes - diag(confusionMatrix)
      if(performanceType == "balanced accuracy")
        mean(diag(confusionMatrix) / classSizes)
      else
        mean(classErrors / classSizes)
    }, actualClasses, predictedClasses, SIMPLIFY = FALSE))
    if(performanceType == "accuracy")
      performanceName <- "Balanced Accuracy"
    else # It is "error".
      performanceName <- "Balanced Error Rate"
  } else { # Metrics for which true positives, true negatives, false positives, false negatives must be calculated.
    performanceName <- switch(performanceType,
                              `micro precision` = "Micro Precision",
                              `micro recall` = "Micro Recall",
                              `micro F1` = "Micro F1 Score",
                              `macro precision` = "Macro Precision",
                              `macro recall` = "Macro Recall",
                              `macro F1` = "Macro F1 Score",
                               matthews = "Matthews Correlation Coefficient")
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
        macroR <- sum(truePositives / (truePositives + falseNegatives)) / nrow(confusionMatrix)
        if(performanceType == "macro recall") return(macroR)
      }
      if(performanceType == "macro F1")
      {
        return(2 * macroP * macroR / (macroP + macroR))
      }
      if(performanceType == "matthews")
      {
        return(unname((truePositives[2] * trueNegatives[2] - falsePositives[2] * falseNegatives[2]) / sqrt((truePositives[2] + falsePositives[2]) * (truePositives[2] + falseNegatives[2]) * (trueNegatives[2] + falsePositives[2]) * (trueNegatives[2] + falseNegatives[2]))))
      }
    }, actualClasses, predictedClasses, SIMPLIFY = FALSE))
  }

  list(name = performanceName, values = performanceValues)
}