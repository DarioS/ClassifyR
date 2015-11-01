setGeneric("calcPerformance", function(result, performanceType, ...)
{standardGeneric("calcPerformance")})

setMethod("calcPerformance", c("ClassifyResult"), function(result, performanceType, ...)
{
  predictions <- lapply(result@predictions, function(sample) sample[, "label"])
  correctClasses <- lapply(result@predictions, function(sample)
    factor(actualClasses(result)[sample[, "sample"]], levels = levels(actualClasses(result)), ordered = TRUE))
  classData <- ROCR::prediction(lapply(predictions, as.numeric), correctClasses)
  if(performanceType == "balanced") # ROCR doesn't do this, currently.
  {
    falseNegativeRate <- sapply(ROCR::performance(classData, "fnr")@y.values, "[[", 2)
    falsePositiveRate <- sapply(ROCR::performance(classData, "fpr")@y.values, "[[", 2)
    falseNegativeRate[is.nan(falseNegativeRate)] <- 0 # When all of the classes are of the second level.
    falsePositiveRate[is.nan(falsePositiveRate)] <- 0 # When all of the classes are of the first level.
    performanceValues <- rowMeans(matrix(c(falseNegativeRate, falsePositiveRate), ncol = 2))
    performanceName <- "Balanced Error Rate"
  } else {
    performanceData <- ROCR::performance(classData, performanceType, ...)
    performanceValues <- sapply(performanceData@y.values, "[[", 2)
    performanceName <- performanceData@y.name
  }
  result@performance[[performanceName]] <- performanceValues
  result
})