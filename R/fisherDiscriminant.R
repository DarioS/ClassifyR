setGeneric("fisherDiscriminant", function(expression, ...)
           {standardGeneric("fisherDiscriminant")})

setMethod("fisherDiscriminant", "matrix", 
          function(expression, classes, ...)
{ 
  colnames(expression) <- NULL # Might be duplicates because of sampling with replacement.            
  features <- rownames(expression)
  groupsTable <- data.frame(class = classes)
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features
  fisherDiscriminant(exprSet, ...)
})

setMethod("fisherDiscriminant", "ExpressionSet", 
          function(expression, test, returnType = c("label", "score", "both"), verbose = 3)
{
  returnType <- match.arg(returnType)
  classes <- pData(expression)[, "class"]
  expression <- exprs(expression)      
  oneClassTraining <- which(classes == levels(classes)[1])
  otherClassTraining <- which(classes == levels(classes)[2])
  varOneClass <- apply(expression[, oneClassTraining], 1, var)
  varOtherClass <- apply(expression[, otherClassTraining], 1, var)
  varAll <- ((length(varOneClass) - 1) * varOneClass + (length(varOtherClass) - 1)
             * varOtherClass) / (length(oneClassTraining) + length(otherClassTraining) - 2)
  aT <- (rowMeans(expression[, oneClassTraining]) - rowMeans(expression[, otherClassTraining])) / varAll
  criticalValue <- 0.5 * aT %*% as.matrix(rowMeans(expression[, oneClassTraining]) +
                                          rowMeans(expression[, otherClassTraining]))
  
  if(verbose == 3)
    message("Critical value calculated.")
  
  labels <- factor(apply(test, 2, function(testSample)
  {
    if(aT %*% as.matrix(testSample) >= criticalValue)
      levels(classes)[1]
    else
      levels(classes)[2]
  }), levels = levels(classes))
  scores <- apply(test, 2, function(testSample) -1 * (aT %*% as.matrix(testSample))) # In reference to the second level of 'classes'. 
  
  switch(returnType, label = labels, score = scores, both = data.frame(label = labels, score = scores))
})