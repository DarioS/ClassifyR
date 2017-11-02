setGeneric("calcPerformance", function(result, performanceType, ...)
{standardGeneric("calcPerformance")})

setMethod("calcPerformance", c("ClassifyResult"), function(result, performanceType, ...)
{
  classLevels <- levels(actualClasses(result))
  predictions <- lapply(result@predictions, function(sample) sample[, "label"])
  correctClasses <- lapply(result@predictions, function(sample)
    factor(actualClasses(result)[sample[, "sample"]], levels = classLevels, ordered = TRUE))
  classData <- ROCR::prediction(lapply(predictions, as.numeric), correctClasses)
  if(performanceType == "balanced") # ROCR can't do this.
  { # Can calculate with two or more class levels.
    performanceValues <- unlist(mapply(function(iterationCorrectClasses, iterationPredictions)
    {
      confusionMatrix <- table(iterationCorrectClasses, iterationPredictions)
      classSizes <- rowSums(confusionMatrix)
      classErrors <- classSizes - diag(confusionMatrix)
      mean(classErrors / classSizes)
    }, correctClasses, predictions, SIMPLIFY = FALSE))
    performanceName <- "Balanced Error Rate"
  } else if(performanceType %in% c("sample error", "sample accuracy"))
  {
    resultTable <- data.frame(do.call(rbind, result@predictions))
    resultTable[, "actual"] <- result@actualClasses[resultTable[, "sample"]]
    allIDs <- 1:length(sampleNames(result))
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
  } else { # Use ROCR calculation. Classes must have exactly two levels.
    performanceData <- ROCR::performance(classData, performanceType, ...)
    performanceValues <- sapply(performanceData@y.values, "[[", 2)
    performanceName <- performanceData@y.name
  }
  result@performance[[performanceName]] <- performanceValues
  result
})