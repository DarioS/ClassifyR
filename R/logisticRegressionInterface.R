setGeneric("logisticRegressionTrainInterface", function(measurements, ...)
{standardGeneric("logisticRegressionTrainInterface")})

setMethod("logisticRegressionTrainInterface", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, ...)
{
  logisticRegressionTrainInterface(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

# Clinical data only.
setMethod("logisticRegressionTrainInterface", "DataFrame", function(measurements, classes, ..., verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)

  if(!requireNamespace("mlogit", quietly = TRUE))
    stop("The package 'mlogit' could not be found. Please install it.")
  if(!requireNamespace("mnlogit", quietly = TRUE))
    stop("The package 'mnlogit' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting logistic regression classifier to data.")

  # Reshape data for input into the fitting function.
  measuredVariables <- colnames(measurements)
  measurements[, "class"] <- classes
  reshaped <- mlogit::mlogit.data(as.data.frame(measurements), choice = "class", shape = "wide")
  modelFormula <- formula(paste("class ~ 1 |", paste(measuredVariables, collapse = '+'), "| 1"))
  
  mnlogit::mnlogit(modelFormula, data = reshaped, ...)
})

# One or more omics datasets, possibly with clinical data.
setMethod("logisticRegressionTrainInterface", "MultiAssayExperiment",
function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, restrict = NULL)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  
  logisticRegressionTrainInterface(measurements, classes, ...)
})

setGeneric("logisticRegressionPredictInterface", function(model, test, ...)
{standardGeneric("logisticRegressionPredictInterface")})

# Matrix of numeric measurements.
setMethod("logisticRegressionPredictInterface", c("mnlogit", "matrix"),
          function(model, test, ...)
{
  logisticRegressionPredictInterface(model, DataFrame(t(test), check.names = FALSE), ...)
})

# Clinical data only.
setMethod("logisticRegressionPredictInterface", c("mnlogit", "DataFrame"), function(model, test, classes = NULL, verbose = 3)
{
  if(!requireNamespace("mlogit", quietly = TRUE))
    stop("The package 'mlogit' could not be found. Please install it.")
  if(!requireNamespace("mnlogit", quietly = TRUE))
    stop("The package 'mnlogit' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained logistic regression classifier.")
  
  if(!is.null(classes)) # Remove them.
  {
    splitDataset <- .splitDataAndClasses(test, classes) # Remove classes, if present.
    test <- splitDataset[["measurements"]]
  }

  # Reshape data for input into the predicting function.
  measuredVariables <- colnames(test)
  test[, "placeholder"] <- model[["choices"]][1] # A hack.
  reshaped <- mlogit::mlogit.data(as.data.frame(test), choice = "placeholder", shape = "wide")
  
  unname(predict(model, reshaped, probability = FALSE))
})

# One or more omics datasets, possibly with clinical data.
setMethod("logisticRegressionPredictInterface", c("mnlogit", "MultiAssayExperiment"),
          function(model, test, targets = names(test), ...)
{
  tablesAndClasses <- .MAEtoWideTable(test, targets, restrict = NULL)
  logisticRegressionPredictInterface(model, tablesAndClasses[["dataTable"]], ...)
})