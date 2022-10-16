# An Interface for ranger Package's randomForest Function
randomForestTrainInterface <- function(measurementsTrain, outcomeTrain, mTryProportion = 0.5, ..., verbose = 3)
{
  if(!requireNamespace("ranger", quietly = TRUE))
    stop("The package 'ranger' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting random forest classifier to training data.")
  mtry <- round(mTryProportion * ncol(measurementsTrain)) # Number of features to try.
      
  # Convert to base data.frame as randomForest doesn't understand DataFrame.
  ranger::ranger(x = as(measurementsTrain, "data.frame"), y = outcomeTrain, mtry = mtry, importance = "impurity_corrected", ...)
}
attr(randomForestTrainInterface, "name") <- "randomForestTrainInterface"
    
# forest is of class ranger
randomForestPredictInterface <- function(forest, measurementsTest, ..., returnType = c("both", "class", "score"), verbose = 3)
{
  returnType <- match.arg(returnType)
  classes <- forest$forest$levels
  if(verbose == 3)
    message("Predicting using random forest.")  
  measurementsTest <- as.data.frame(measurementsTest)
  classPredictions <- predict(forest, measurementsTest)$predictions
  classScores <- predict(forest, measurementsTest, predict.all = TRUE)[[1]]
  classScores <- t(apply(classScores, 1, function(sampleRow) table(factor(classes[sampleRow], levels = classes)) / forest$forest$num.trees))
  switch(returnType, class = classPredictions,
         score = classScores,
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
}

################################################################################
#
# Get selected features
#
################################################################################

forestFeatures <- function(forest)
                  {
                    rankedFeaturesIndices <- order(ranger::importance(forest), decreasing = TRUE)
                    selectedFeaturesIndices <- which(ranger::importance(forest) > 0)
                    list(rankedFeaturesIndices, selectedFeaturesIndices)
                  }