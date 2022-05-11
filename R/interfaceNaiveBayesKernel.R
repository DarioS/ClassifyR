#' Classification Using A Bayes Classifier with Kernel Density Estimates
#' 
#' Kernel density estimates are fitted to the training data and a naive Bayes
#' classifier is used to classify samples in the test data.
#' 
#' If \code{difference} is \code{"weighted"}, then a sample's predicted class
#' is the class with the largest sum of weights, each scaled for the number of
#' samples in the training data of each class. Otherwise, when
#' \code{difference} is \code{"unweighted"}, each feature has an equal vote,
#' and votes for the class with the largest weight, scaled for class sizes in
#' the training set.
#' 
#' The variable name of each feature's measurements in the iteration over all
#' features is \code{featureValues}.  This is important to know if each
#' feature's measurements need to be referred to in the specification of
#' \code{densityParameters}, such as for specifying the range of x values of
#' the density function to be computed.  For example, see the default value of
#' \code{densityParameters} above.
#' 
#' If \code{weight} is \code{"crossover distance"}, the crossover points are
#' computed by considering the distance between y values of all of the
#' densities at every x value. x values for which a class density crosses any
#' other class' density are used as the crossover points for that class.
#' 
#' @aliases naiveBayesKernel naiveBayesKernel,matrix-method
#' naiveBayesKernel,DataFrame-method
#' naiveBayesKernel,MultiAssayExperiment-method
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
#' If of type \code{\link{DataFrame}} or \code{\link{MultiAssayExperiment}}, the data set is subset
#' to only those features of type \code{numeric}.
#' @param classesTrain A vector of class labels of class \code{\link{factor}} of the
#' same length as the number of samples in \code{measurementsTrain} if it is a
#' \code{\link{matrix}} or a \code{\link{DataFrame}} or a character vector of length 1
#' containing the column name in \code{measurementsTrain} if it is a \code{\link{DataFrame}} or the
#' column name in \code{colData(measurementsTrain)} if \code{measurementsTrain} is a
#' \code{\link{MultiAssayExperiment}}. If a column name, that column will be
#' removed before training.
#' @param measurementsTest An object of the same class as \code{measurementsTrain} with no
#' samples in common with \code{measurementsTrain} and the same number of features
#' as it.
#' @param targets If \code{measurementsTrain} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"sampleInfo"} is also a valid value
#' and specifies that numeric variables from the sample information data table will be
#' used.
#' @param ... Unused variables by the three top-level methods passed to the
#' internal method which does the classification.
#' @param densityFunction Default: \code{\link{density}}. A function which will
#' return a probability density, which is essentially a list with x and y
#' coordinates.
#' @param densityParameters A list of options for \code{densityFunction}.
#' Default: \code{list(bw = "nrd0", n = 1024, from =
#' expression(min(featureValues)), to = expression(max(featureValues))}.
#' @param difference Default: \code{"unweighted"}. Either \code{"unweighted"},
#' \code{"weighted"}. In weighted mode, the difference in densities is summed
#' over all features. If unweighted mode, each feature's vote is worth the
#' same.
#' @param weighting Default: \code{"height difference"}. Either \code{"height difference"}
#' or \code{"height difference"}. The type of weight to calculate. For
#' \code{"height difference"}, the weight of each prediction is equal to the
#' vertical distance between the highest density and the second-highest, for a
#' particular value of x. For \code{"crossover distance"}, the x positions
#' where two densities cross is firstly calculated.  The predicted class is the
#' class with the highest density at the particular value of x and the weight
#' is the distance of x from the nearest density crossover point.
#' @param minDifference Default: 0. The minimum difference in density height
#' between the highest density and second-highest for a feature to be allowed
#' to vote. If no features for a particular sample have a difference large enough,
#' the class predicted is simply the largest class.
#' @param returnType Default: \code{"both"}. Either \code{"class"},
#' \code{"score"} or \code{"both"}.  Sets the return value from the prediction
#' to either a vector of predicted classes, a matrix of scores with columns
#' corresponding to classes, as determined by the factor levels of
#' \code{classes}, or both a column of predicted classes and columns of class
#' scores in a \code{data.frame}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return A vector or data frame of class prediction information (i.e. classes
#' and/or scores), as long as the number of samples in the test data.
#' @author Dario Strbenac, John Ormerod
#' @examples
#' 
#'   trainMatrix <- matrix(rnorm(1000, 8, 2), nrow = 10)
#'   classesTrain <- factor(rep(c("Poor", "Good"), each = 5))
#'   rownames(trainMatrix) <- paste("Sample", 1:10)
#'   
#'   # Make first 30 genes increased in value for poor samples.
#'   trainMatrix[1:5, 1:30] <- trainMatrix[1:5, 1:30] + 5
#'   
#'   testMatrix <- matrix(rnorm(1000, 8, 2), nrow = 10)
#'   rownames(testMatrix) <- paste("Sample", 11:20)
#'   
#'   # Make first 30 genes increased in value for sixth to tenth samples.
#'   testMatrix[6:10, 1:30] <- testMatrix[6:10, 1:30] + 5
#'   
#'   naiveBayesKernel(trainMatrix, classesTrain, testMatrix)
#' 
#' @rdname naiveBayesKernel
#' @usage NULL
#' @export
setGeneric("naiveBayesKernel", function(measurementsTrain, ...)
           standardGeneric("naiveBayesKernel"))

#' @rdname naiveBayesKernel
#' @export
setMethod("naiveBayesKernel", "matrix", # Matrix of numeric measurements.
          function(measurementsTrain, classesTrain, measurementsTest, ...)
{
  naiveBayesKernel(S4Vectors::DataFrame(measurementsTrain, check.names = FALSE),
                   classesTrain,
                   S4Vectors::DataFrame(measurementsTest, check.names = FALSE), ...)
})

#' @rdname naiveBayesKernel
#' @export
setMethod("naiveBayesKernel", "DataFrame", # Sample information data or one of the other inputs, transformed.
          function(measurementsTrain, classesTrain, measurementsTest,
                   densityFunction = density, densityParameters = list(bw = "nrd0", n = 1024, from = expression(min(featureValues)), to = expression(max(featureValues))),
                   difference = c("unweighted", "weighted"),
                   weighting = c("height difference", "crossover distance"),
                   minDifference = 0, returnType = c("both", "class", "score"), verbose = 3)
{
  splitDataset <- .splitDataAndOutcomes(measurementsTrain, classesTrain)
  trainingMatrix <- splitDataset[["measurements"]]
  classesTrain <- splitDataset[["outcomes"]]
  testingMatrix <- as.matrix(measurementsTest[, colnames(trainingMatrix), drop = FALSE])
  
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  
  difference <- match.arg(difference)
  weighting <- match.arg(weighting)
  returnType <- match.arg(returnType)
  
  classesSizes <- sapply(levels(classesTrain), function(class) sum(classesTrain == class))
  largestClass <- names(classesSizes)[which.max(classesSizes)[1]]
  
  if(verbose == 3)
    message("Fitting densities.")
  
  featuresDensities <- lapply(measurementsTrain, function(featureValues)
  {
    densityParameters <- lapply(densityParameters, function(parameter) eval(parameter))
    lapply(levels(classesTrain), function(class)
    {
      aClassMeasurements <- featureValues[classesTrain == class]  
      do.call(densityFunction, c(list(aClassMeasurements), densityParameters))
    }) # A fitted density for each class.
  })

  classesScaleFactors <- classesSizes / nrow(trainingMatrix)
  splines <- lapply(featuresDensities, function(featureDensities) 
             {
               mapply(function(featureDensity, scaleFactor)
               {
                 splinefun(featureDensity[['x']], featureDensity[['y']] * scaleFactor, "natural")
               }, featureDensities, classesScaleFactors)
             })
  
  if(verbose == 3)
    message("Calculating vertical distances between class densities.")

  # Needed even if horizontal distance weighting is used to determine the predicted class.
  posteriorsVertical <- mapply(function(featureSplines, testSamples)
  {
    sapply(1:length(levels(classesTrain)), function(classIndex)
    {
      featureSplines[[classIndex]](testSamples)
    })
  }, splines, as.data.frame(testingMatrix), SIMPLIFY = FALSE)
    
  classesVertical <- sapply(posteriorsVertical, function(featureVertical)
  {
      apply(featureVertical, 1, function(sampleVertical) levels(classesTrain)[which.max(sampleVertical)])
  }) # Matrix, rows are test samples, columns are features.
    
  distancesVertical <- sapply(posteriorsVertical, function(featureVertical)
  { # Vertical distance between highest density and second-highest, at a particular value.
    apply(featureVertical, 1, function(sampleVertical)
    {
      twoHighest <- sort(sampleVertical, decreasing = TRUE)[1:2]
      Reduce('-', twoHighest)
    })
  }) # Matrix, rows are test samples, columns are features.
  
  if(difference == "weighted" && weight == "crossover distance")
  {
    if(verbose == 3)
      message("Calculating horizontal distances to crossover points of class densities.")
 
    classesVerticalIndices <- matrix(match(classesVertical, levels(classesTrain)),
                                     nrow = nrow(classesVertical), ncol = ncol(classesVertical))
    distancesHorizontal <- mapply(function(featureDensities, testSamples, predictedClasses)
    {
      classesCrosses <- .densitiesCrossover(featureDensities)
      classesDistances <- sapply(classesCrosses, function(classCrosses)
      {
        sapply(testSamples, function(testSample) min(abs(testSample - classCrosses)))
      })
      classesDistances[cbind(1:nrow(classesDistances), predictedClasses)]
    }, featuresDensities, test, as.data.frame(classesVerticalIndices)) # Matrix of horizontal distances to nearest cross-over involving the predicted class.
  }

  if(verbose == 3)
  {
    switch(returnType, class = message("Determining class labels."),
                       both = message("Calculating class scores and determining class labels."),
                       score = message("Calculating class scores."))
  }
  
  allDistances <- switch(weighting, `height difference` = distancesVertical,
                                    `crossover distance` = distancesHorizontal)

  predictions <- do.call(rbind, lapply(1:nrow(allDistances), function(sampleRow)
  {
    useFeatures <- abs(allDistances[sampleRow, ]) > minDifference
    if(all(useFeatures == FALSE)) # No features have a large enough density difference.
    {                          # Simply vote for the larger class.
      classPredicted <- largestClass
      classScores <- classesSizes / length(classesTrain)
    } else { # One or more features are available to vote with.
      distancesUsed <- allDistances[sampleRow, useFeatures]
      classPredictionsUsed <- factor(classesVertical[sampleRow, useFeatures], levels(classesTrain))
      if(difference == "unweighted")
      {
        classScores <- table(classPredictionsUsed)
        classScores <- setNames(as.vector(classScores), levels(classesTrain))
      } else { # Weighted voting.
        classScores <- tapply(distancesUsed, classPredictionsUsed, sum)
        classScores[is.na(classScores)] <- 0
      }
      classScores <- classScores / sum(classScores) # Make different feature selection sizes comparable.
      classPredicted <- names(classScores)[which.max(classScores)]
    }

    data.frame(class = factor(classPredicted, levels = levels(classesTrain)), t(classScores), check.names = FALSE)
  }))

  switch(returnType, class = predictions[, "class"],
         score = predictions[, 2:ncol(predictions)],
         both = predictions)
})

#' @rdname naiveBayesKernel
#' @export
setMethod("naiveBayesKernel", "MultiAssayExperiment",
          function(measurementsTrain, measurementsTest, targets = names(measurements), classesTrain, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurementsTrain, targets, classesTrain)
  trainingMatrix <- tablesAndClasses[["dataTable"]]
  classesTrain <- tablesAndClasses[["outcomes"]]
  testingMatrix <- .MAEtoWideTable(measurementsTest, targets)
            
  .checkVariablesAndSame(trainingMatrix, testingMatrix)
  naiveBayesKernel(trainingMatrix, classesTrain, testingMatrix, ...)
})