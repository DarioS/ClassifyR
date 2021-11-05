setGeneric("SVMtrainInterface", function(measurements, ...)
standardGeneric("SVMtrainInterface"))

setMethod("SVMtrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  SVMtrainInterface(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

# Clinical data or one of the other inputs, transformed.
setMethod("SVMtrainInterface", "DataFrame", function(measurements, classes, ..., verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  trainingMatrix <- as.matrix(splitDataset[["measurements"]])
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]

  if(!requireNamespace("e1071", quietly = TRUE))
    stop("The package 'e1071' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting SVM classifier to data.")

  e1071::svm(measurements, classes, probability = TRUE, ...)
})

setMethod("SVMtrainInterface", "MultiAssayExperiment",
function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  
  if(ncol(measurements) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    SVMtrainInterface(measurements, classes, ...)
})

setGeneric("SVMpredictInterface", function(model, test, ...)
standardGeneric("SVMpredictInterface"))

setMethod("SVMpredictInterface", c("svm", "matrix"),
          function(model, test, ...)
{
  SVMpredictInterface(model, DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("SVMpredictInterface", c("svm", "DataFrame"), function(model, test, classes = NULL, returnType = c("both", "class", "score"), classifierName = "Support Vector Machine", verbose = 3)
{
  returnType <- match.arg(returnType)
  if(!is.null(classes))
  {
    splitDataset <- .splitDataAndClasses(test, classes) # Remove classes, if present.
    testMatrix <- splitDataset[["measurements"]]
  } else {testMatrix <- test}
  isNumeric <- sapply(testMatrix, is.numeric)
  testMatrix <- testMatrix[, isNumeric, drop = FALSE]
  
  if(!requireNamespace("e1071", quietly = TRUE))
    stop("The package 'e1071' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained SVM classifier.")
  
  classPredictions <- predict(model, test, probability = TRUE)
  classScores <- attr(classPredictions, "probabilities")[, model[["levels"]], drop = FALSE]
  attr(classPredictions, "probabilities") <- NULL
  switch(returnType, class = classPredictions, score = classScores,
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
})

setMethod("SVMpredictInterface", c("svm", "MultiAssayExperiment"),
          function(model, test, targets = names(test), ...)
{
  tablesAndClasses <- .MAEtoWideTable(test, targets)
  test <- tablesAndClasses[["dataTable"]] # Remove any classes, if present.
            
  if(ncol(test) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    SVMpredictInterface(model, test, ...)
})