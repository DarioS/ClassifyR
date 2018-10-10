setGeneric("classifyInterface", function(measurements, ...)
{standardGeneric("classifyInterface")})

setMethod("classifyInterface", "matrix", # Matrix of integer measurements.
          function(measurements, classes, test, ...)
{
  classifyInterface(DataFrame(t(measurements), check.names = FALSE),
                    classes,
                    DataFrame(t(test), check.names = FALSE), ...)
})

# Clinical data or one of the other inputs, transformed.
setMethod("classifyInterface", "DataFrame", function(measurements, classes, test, ..., returnType = c("class", "score", "both"), verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  trainingMatrix <- splitDataset[["measurements"]]
  isInteger <- apply(trainingMatrix, 2, is.integer)
  trainingMatrix <- as.matrix(trainingMatrix[, isInteger, drop = FALSE])
  isInteger <- sapply(test, is.integer)
  testingMatrix <- as.matrix(test[, isInteger, drop = FALSE])
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  
  returnType <- match.arg(returnType)
  
  if(!requireNamespace("PoiClaClu", quietly = TRUE))
    stop("The package 'PoiClaClu' could not be found. Please install it.")
  
  if(verbose == 3)
    message("Fitting Poisson LDA classifier to training data and making predictions on test
            data.")

  predicted <- PoiClaClu::Classify(trainingMatrix, classes, testingMatrix, ...)
  classPredictions <- predicted[["ytehat"]]
  classScores <- predicted[["discriminant"]]
  colnames(classScores) <- levels(classes)
  switch(returnType, class = classPredictions, # Factor vector.
         score = classScores, # Numeric matrix.
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
})

setMethod("classifyInterface", "MultiAssayExperiment",
function(measurements, test, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, "integer")
  trainingMatrix <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  testingMatrix <- .MAEtoWideTable(test, targets, "integer")
            
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  classifyInterface(trainingMatrix, classes, testingMatrix, ...)
})