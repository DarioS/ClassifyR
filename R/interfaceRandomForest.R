# An Interface for randomForest Package's randomForest Function
randomForestTrainInterface <- function(measurementsTrain, classesTrain, mTryProportion = 0.5, ..., verbose = 3)
{
  if(!requireNamespace("randomForest", quietly = TRUE))
    stop("The package 'randomForest' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting random forest classifier to training data and making predictions on test
            data.")
  mtry <- round(mTryProportion * ncol(measurementsTrain)) # Number of features to try.
      
  # Convert to base data.frame as randomForest doesn't understand DataFrame.
  randomForest::randomForest(as(measurementsTrain, "data.frame"), classesTrain, mtry = mtry, keep.forest = TRUE, ...)
}
attr(randomForestTrainInterface, "name") <- "randomForestTrainInterface"
    
# forest is of class randomForest
randomForestPredictInterface <- function(forest, measurementsTest, ..., returnType = c("both", "class", "score"), verbose = 3)
{
  returnType <- match.arg(returnType)
  if(verbose == 3)
    message("Predicting using random forest.")  
  measurementsTest <- as.data.frame(measurementsTest)
  classPredictions <- predict(forest, measurementsTest)
  classScores <- predict(forest, measurementsTest, type = "vote")[, forest[["classes"]], drop = FALSE]
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
                    rankedFeaturesIndices <- order(randomForest::importance(forest), decreasing = TRUE)
                    selectedFeaturesIndices <- randomForest::varUsed(forest, count = FALSE)
                    list(rankedFeaturesIndices, selectedFeaturesIndices)
                  }