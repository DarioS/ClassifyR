# An Interface for glmnet Package's glmnet Function. Generalised linear models with sparsity.

elasticNetGLMtrainInterface <- function(measurementsTrain, classesTrain, lambda = NULL, ..., verbose = 3)
{
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting elastic net regularised GLM classifier to data.")
    
  measurementsMatrix <- glmnet::makeX(as(measurementsTrain, "data.frame"))

  fitted <- glmnet::glmnet(measurementsMatrix, classesTrain, family = "multinomial", ...)

  if(is.null(lambda)) # fitted has numerous models for automatically chosen lambda values.
  { # Pick one lambda based on resubstitution performance.
    bestLambda <- fitted[["lambda"]][which.min(sapply(fitted[["lambda"]], function(lambda) # Largest Lambda with minimum balanced error rate.
    {
      classPredictions <- factor(as.character(predict(fitted, measurementsMatrix, s = lambda, type = "class")), levels = fitted[["classnames"]])
      calcExternalPerformance(classesTrain, classPredictions, "Balanced Error")
    }))[1]]
    attr(fitted, "tune") <- list(lambda = bestLambda)
  }
  fitted
}
attr(elasticNetGLMtrainInterface, "name") <- "elasticNetGLMtrainInterface"

# model is of class multnet
elasticNetGLMpredictInterface <- function(model, measurementsTest, lambda, ..., returnType = c("both", "class", "score"), verbose = 3)
{ # ... just consumes emitted tuning variables from .doTrain which are unused.
  returnType <- match.arg(returnType)
  
  # Ensure that testing data has same columns names in same order as training data.
  # Remove those annoying backquotes which glmnet adds if variables have spaces in names.
  measurementsTest <- measurementsTest[, gsub('`', '', rownames(model[["beta"]][[1]]))]
  
  if(!requireNamespace("glmnet", quietly = TRUE))
    stop("The package 'glmnet' could not be found. Please install it.")
  if(verbose == 3)
    message("Predicting classes using trained elastic net regularised GLM classifier.")

  if(missing(lambda)) # Tuning parameters are not passed to prediction functions.
    lambda <- attr(model, "tune")[["lambda"]] # Sneak it in as an attribute on the model.

  testMatrix <- glmnet::makeX(as(measurementsTest, "data.frame"))
  testMatrix <- testMatrix[, rownames(model[["beta"]][[1]])]
  
  classPredictions <- factor(as.character(predict(model, testMatrix, s = lambda, type = "class")), levels = model[["classnames"]])
  classScores <- predict(model, testMatrix, s = lambda, type = "response")[, , 1]
  
  if(is.matrix(classScores))
    classScores <- classScores[, model[["classnames"]]]
  else # Leave-one-out cross-validation likely used and glmnet doesn't have consistent return types.
    classScores <- t(classScores[model[["classnames"]]])
  
  switch(returnType, class = classPredictions, # Factor vector.
         score = classScores, # Numeric matrix.
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
}

################################################################################
#
# Get selected features (i.e. non-zero model coefficients)
#
################################################################################

elasticNetFeatures <- function(model)
                      {
                        inputFeatures <- rownames(model[["beta"]][[1]])            
                        # Floating point numbers test for equality.
                        whichCoefficientColumn <- which(abs(model[["lambda"]] - attr(model, "tune")[["lambda"]]) < 0.00001)[1]
                        coefficientsUsed <- sapply(model[["beta"]], function(classCoefficients) classCoefficients[, whichCoefficientColumn])
                        featureScores <- rowSums(abs(coefficientsUsed))
                        rankedFeaturesIndices <- order(featureScores, decreasing = TRUE)
                        selectedFeaturesIndices <- which(featureScores != 0)
                        list(rankedFeaturesIndices, selectedFeaturesIndices)
                      }