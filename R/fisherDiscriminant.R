setGeneric("fisherDiscriminant", function(measurements, ...)
           standardGeneric("fisherDiscriminant"))

setMethod("fisherDiscriminant", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, test, ...)
{
  fisherDiscriminant(DataFrame(t(measurements[, , drop = FALSE]), check.names = FALSE),
                     classes,
                     DataFrame(t(test[, , drop = FALSE]), check.names = FALSE), ...)
})

setMethod("fisherDiscriminant", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, test, returnType = c("both", "class", "score"), verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  trainingMatrix <- splitDataset[["measurements"]]
  isNumeric <- sapply(trainingMatrix, is.numeric)
  trainingMatrix <- as.matrix(trainingMatrix[, isNumeric, drop = FALSE])
  isNumeric <- sapply(test, is.numeric)
  testingMatrix <- as.matrix(test[, isNumeric, drop = FALSE])
            
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  returnType <- match.arg(returnType)

  oneClassTraining <- which(classes == levels(classes)[1])
  otherClassTraining <- which(classes == levels(classes)[2])
  varOneClass <- apply(measurements[oneClassTraining, ], 2, var)
  varOtherClass <- apply(measurements[otherClassTraining, ], 2, var)
  varAll <- ((length(varOneClass) - 1) * varOneClass + (length(varOtherClass) - 1)
             * varOtherClass) / (length(oneClassTraining) + length(otherClassTraining) - 2)
  aT <- (apply(measurements[oneClassTraining, ], 2, mean) - apply(measurements[otherClassTraining, ], 2, mean)) / varAll
  criticalValue <- 0.5 * aT %*% as.matrix(apply(measurements[oneClassTraining, ], 2, mean) +
                                          apply(measurements[otherClassTraining, ], 2, mean)
                                         )
  
  if(verbose == 3)
    message("Critical value calculated.")
  
  classes <- factor(apply(test, 1, function(testSample)
  {
    if(aT %*% as.matrix(testSample) >= criticalValue)
      levels(classes)[1]
    else
      levels(classes)[2]
  }), levels = levels(classes))
  scores <- apply(test, 1, function(testSample) -1 * (aT %*% as.matrix(testSample))) # In reference to the second level of 'classes'. 
  
  switch(returnType, class = classes,
                     score = scores,
                     both = data.frame(class = classes, score = scores))  
})

# One or more omics data sets, possibly with clinical data.
setMethod("fisherDiscriminant", "MultiAssayExperiment", 
          function(measurements, test, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  trainingMatrix <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  testingMatrix <- .MAEtoWideTable(test, targets)
  
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  fisherDiscriminant(trainingMatrix, classes, testingMatrix, ...)
})