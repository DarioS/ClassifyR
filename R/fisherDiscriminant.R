setGeneric("fisherDiscriminant", function(measurements, ...)
           {standardGeneric("fisherDiscriminant")})

setMethod("fisherDiscriminant", "matrix", 
          function(measurements, classes, test, ...)
{
  .fisherDiscriminant(DataFrame(t(measurements[, , drop = FALSE]), check.names = FALSE), classes,
                      DataFrame(t(test[, , drop = FALSE]), check.names = FALSE), ...)
})

setMethod("fisherDiscriminant", "DataFrame", 
          function(measurements, classes, test, ...)
          {
            if(ncol(measurements) != ncol(test))
              stop("Training data table and testing data table have a different number of features.")
            splitDataset <- .splitDataAndClasses(measurements, classes)
            measurements <- splitDataset[["measurements"]]
            isNumeric <- apply(measurements, 2, is.numeric)
            measurements <- measurements[, isNumeric, drop = FALSE]
            test <- test[, isNumeric, drop = FALSE]
            if(sum(isNumeric) == 0)
              stop("All features are not numeric but at least one must be.")
            .fisherDiscriminant(measurements, splitDataset[["classes"]], test, ...)
          })

setMethod("fisherDiscriminant", "MultiAssayExperiment", 
          function(measurements, test, targets = names(measurements))
{
  if(!all(targets %in% c(names(measurements), "colData")))
    stop("Some table names in 'targets' are not assay names in 'measurements' or \"colData\".") 
            
  measurements <- measurements[, , targets]
  trainingMatrix <- wideFormat(measurements, colDataCols = colnames(colData(measurements)),
                               check.names = FALSE)
  classes <- trainingMatrix[, "class"]
  trainingMatrix <- trainingMatrix[, -match(c("primary", "class"), colnames(trainingMatrix))]
  isNumeric <- apply(trainingMatrix, 2, trainingMatrix)
  trainingMatrix <- trainingMatrix[, isNumeric, drop = FALSE]
  
  testingMatrix <- wideFormat(test, colDataCols = colnames(colData(measurements)),
                              check.names = FALSE)
  primaryKey <- match("primary", colnames(testingMatrix))
  rownames(testingMatrix) <- testingMatrix[, primaryKey]
  testingMatrix <- testingMatrix[, -primaryKey]
  isNumeric <- apply(testingMatrix, 2, testingMatrix)
  testingMatrix <- testingMatrix[, isNumeric, drop = FALSE]
  
  if(ncol(trainingMatrix) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else if(ncol(trainingMatrix) != ncol(testingMatrix))
    stop("Training dataset and testing dataset contain differing numbers of features.")
  else
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