setGeneric("fisherDiscriminant", function(measurements, ...)
           {standardGeneric("fisherDiscriminant")})

setMethod("fisherDiscriminant", "matrix", 
          function(measurements, classes, test, ...)
{
  .fisherDiscriminant(DataFrame(t(measurements[, , drop = FALSE]), check.names = FALSE),
                      classes,
                      DataFrame(t(test[, , drop = FALSE]), check.names = FALSE), ...)
})

setMethod("fisherDiscriminant", "DataFrame", 
          function(measurements, classes, test, ...)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  trainingMatrix <- as.matrix(splitDataset[["measurements"]])
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  isNumeric <- sapply(test, is.numeric)
  testingMatrix <- as.matrix(test[, isNumeric, drop = FALSE])
            
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  .fisherDiscriminant(trainingMatrix, splitDataset[["classes"]], testingMatrix, ...)
})

setMethod("fisherDiscriminant", "MultiAssayExperiment", 
          function(measurements, test, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  trainingMatrix <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  testingMatrix <- .MAEtoWideTable(test, targets)
  
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  .fisherDiscriminant(trainingMatrix, classes, testingMatrix, ...)
})

.fisherDiscriminant <- function(measurements, classes, test, returnType = c("label", "score", "both"), verbose = 3)
{
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
  
  labels <- factor(apply(test, 1, function(testSample)
  {
    if(aT %*% as.matrix(testSample) >= criticalValue)
      levels(classes)[1]
    else
      levels(classes)[2]
  }), levels = levels(classes))
  scores <- apply(test, 1, function(testSample) -1 * (aT %*% as.matrix(testSample))) # In reference to the second level of 'classes'. 
  
  switch(returnType, label = labels, score = scores, both = data.frame(label = labels, score = scores))
}