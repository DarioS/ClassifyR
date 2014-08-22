setGeneric("calcError", function(result, errorType, ...)
{standardGeneric("calcError")})

setMethod("calcError", c("ClassifyResult"), function(result, errorType, ...)
{
  predictions <- lapply(result@predictions, function(sample) sample[, "predicted"])
  correctClasses <- lapply(result@predictions, function(sample)
    factor(actualClasses(result)[sample[, "sample"]], ordered = TRUE))
  classData <- prediction(lapply(predictions, as.numeric), correctClasses)
  if(errorType == "balanced") # ROCR doesn't do this, currently.
  {
    faleNegativeRate <- sapply(performance(classData, "fnr")@y.values, "[[", 2)
    falePositiveRate <- sapply(performance(classData, "fpr")@y.values, "[[", 2)
    errorValues <- rowMeans(matrix(c(faleNegativeRate, falePositiveRate), ncol = 2))
    errorName <- "Balanced Error Rate"
  } else {
    errorData <- performance(classData, errorType, ...)
    errorValues <- sapply(errorData@y.values, "[[", 2)
    errorName <- errorData@y.name
  }
  result@errors[[errorName]] <- errorValues
  result
})