setGeneric("randomForestTrainInterface", function(measurements, ...)
{standardGeneric("randomForestTrainInterface")})
setGeneric("randomForestPredictInterface", function(models, test, ...)
           {standardGeneric("randomForestPredictInterface")})

setMethod("randomForestTrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  randomForestTrainInterface(DataFrame(t(measurements), check.names = FALSE),
                             classes, ...)
})

# Clinical data or one of the other inputs, transformed.
setMethod("randomForestTrainInterface", "DataFrame", function(measurements, classes, ..., verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)

  if(!requireNamespace("randomForest", quietly = TRUE))
    stop("The package 'randomForest' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting random forest classifier to training data and making predictions on test
            data.")

  # Convert to base data.frame as randomForest doesn't understand DataFrame.
  randomForest::randomForest(as.data.frame(splitDataset[["measurements"]]), splitDataset[["classes"]], keep.forest = TRUE, ...)
})

setMethod("randomForestTrainInterface", "MultiAssayExperiment",
function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, restrict = NULL)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  
  randomForestTrainInterface(measurements, classes, ...)
})

setGeneric("randomForestPredictInterface", function(forest, test, ...)
           {standardGeneric("randomForestPredictInterface")})

setMethod("randomForestPredictInterface", c("randomForest", "matrix"), function(forest, test, ...)
{
  randomForestPredictInterface(forest, DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("randomForestPredictInterface", c("randomForest", "DataFrame"),
function(forest, test, ..., returnType = c("class", "score", "both"), verbose = 3)
{
  returnType <- match.arg(returnType)
  if(verbose == 3)
    message("Predicting using random forest.")  
  
  classPredictions <- predict(forest, test)
  classScores <- predict(forest, test, type = "vote")[, forest[["classes"]], drop = FALSE]
  switch(returnType, class = classPredictions,
         score = classScores,
         both = data.frame(class = classPredictions, classScores))
})

# One or more omics data sets, possibly with clinical data.
setMethod("randomForestPredictInterface", c("randomForest", "MultiAssayExperiment"),
          function(forest, test, targets = names(test), ...)
{
  testingTable <- .MAEtoWideTable(test, targets)
  randomForestPredictInterface(models, testingTable, ...)
})