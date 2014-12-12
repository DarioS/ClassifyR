setClass("TransformParams", representation(
  transform = "function",
  intermediate = "character",  
  otherParams = "list")
)

setGeneric("TransformParams", function(transform, ...)
{standardGeneric("TransformParams")})

setMethod("TransformParams", c("function"),
          function(transform, intermediate = character(0), ...)
          {
            new("TransformParams", transform = transform,
                intermediate = intermediate, otherParams = list(...))
          })

setClassUnion("functionOrList", c("function", "list"))

setClass("ResubstituteParams", representation(
  nFeatures = "numeric",
  performanceType = "character",
  better = "character",
  otherParams = "list")
)

setGeneric("ResubstituteParams", function(nFeatures, performanceType, better = c("lower", "higher"), ...)
{standardGeneric("ResubstituteParams")})

setMethod("ResubstituteParams", numeric(0), function()
{
  new("ResubstituteParams", nFeatures = seq(100, 500, 100), performanceType = "balanced",
      better = "lower")
})

setMethod("ResubstituteParams", c("numeric", "character", "character"),
          function(nFeatures, performanceType, better, ...)
          {
            new("ResubstituteParams", nFeatures = nFeatures, performanceType = performanceType,
                better = better, otherParams = list(...))
          })

setClass("SelectionParams", representation(
  featureSelection = "functionOrList",
  minPresence = "numeric",
  intermediate = "character",
  subsetExpressionData = "logical",  
  otherParams = "list")
)

setGeneric("SelectionParams", function(featureSelection, ...)
{standardGeneric("SelectionParams")})
setMethod("SelectionParams", character(0), function()
{
  new("SelectionParams", featureSelection = limmaSelection, minPresence = 1,
      intermediate = character(0), subsetExpressionData = TRUE,
      otherParams = list(resubstituteParams = ResubstituteParams()))
})
setMethod("SelectionParams", c("functionOrList"),
          function(featureSelection, minPresence = 1, intermediate = character(0),
                   subsetExpressionData = TRUE, ...)
          {
            new("SelectionParams", featureSelection = featureSelection,
                minPresence = minPresence, intermediate = intermediate,
                subsetExpressionData = subsetExpressionData, otherParams = list(...))
          })

setClass("TrainParams", representation(
  classifier = "function",
  transposeExpression = "logical",
  doesTests = "logical",
  intermediate = "character",  
  otherParams = "list")
)

setGeneric("TrainParams", function(classifier, ...)
{standardGeneric("TrainParams")})
setMethod("TrainParams", character(0), function()
{
  new("TrainParams", classifier = dlda, transposeExpression = TRUE, intermediate = character(0),
      doesTests = FALSE)
})
setMethod("TrainParams", c("function"),
          function(classifier, transposeExpression, doesTests, intermediate = character(0), ...)
          {
            new("TrainParams", classifier = classifier, transposeExpression = transposeExpression,
                doesTests = doesTests, intermediate = intermediate, otherParams = list(...))
          })

setClass("PredictParams", representation(
  predictor = "function",
  transposeExpression = "logical",
  multipleResults = "logical",
  intermediate = "character",    
  getClasses = "function",
  otherParams = "list")
)

setGeneric("PredictParams", function(predictor, ...)
{standardGeneric("PredictParams")})
setMethod("PredictParams", character(0), function()
{
  new("PredictParams", predictor = predict, transposeExpression = TRUE, multipleResults = FALSE,
      intermediate = character(0), getClasses = function(result){result[["class"]]})
})
setMethod("PredictParams", c("function"),
          function(predictor, transposeExpression, multipleResults, intermediate = character(0), getClasses, ...)
          {
            new("PredictParams", predictor = predictor, transposeExpression = transposeExpression,
                multipleResults = multipleResults, intermediate = intermediate,
                getClasses = getClasses, otherParams = list(...))
          })

setGeneric("ClassifyResult", function(datasetName, classificationName, originalNames, originalFeatures, ...)
{standardGeneric("ClassifyResult")})
setClass("ClassifyResult", representation(
  datasetName = "character",
  classificationName = "character",
  originalNames = "character",
  originalFeatures = "character",
  chosenFeatures = "list",
  rankedFeatures = "list",
  actualClasses = "factor",
  predictions = "list",
  validation = "list",  
  performance = "list")
)
setMethod("ClassifyResult", c("character", "character", "character", "character"),
          function(datasetName, classificationName, originalNames, originalFeatures,
                   rankedFeatures, chosenFeatures, predictions, actualClasses, validation)
          {
            new("ClassifyResult", datasetName = datasetName, classificationName = classificationName,
                predictions = predictions, rankedFeatures = rankedFeatures, chosenFeatures = chosenFeatures,
                actualClasses = actualClasses, validation = validation,
                originalNames = originalNames, originalFeatures = originalFeatures)
          })
setMethod("show", c("ClassifyResult"),
          function(object)
          {
            cat("An object of class 'ClassifyResult'.\n")
            cat("Dataset Name: ", object@datasetName, ".\n", sep = '')
            cat("Classification Name: ", object@classificationName, ".\n", sep = '')
            cat("Validation: ") 
            if(object@validation[[1]] == "leave")
              cat("Leave ", object@validation[[2]], " out cross-validation.\n", sep = '')
            else if(object@validation[[1]] == "split")
              cat("Split cross-validation, ", object@validation[[2]], " percent of samples in test set.",
                  sep = '')
            else
              cat(object@validation[[3]], " fold cross-validation of ", object@validation[[2]],
                  " resamples.\n", sep = '')
            cat("Predictions: List of data frames of length ", length(object@predictions),
                ".\n", sep = '')
            if(object@validation[[1]] %in% c("split", "leave"))
              cat("Features: List of length ", length(object@chosenFeatures), " of row indices.\n", sep = '')
            else
              cat("Features: List of length ", length(object@chosenFeatures), " of lists of length ",
                  length(object@chosenFeatures[[1]]), " of row indices.\n", sep = '')
            if(length(object@performance) > 0)
              cat("Performance Measures: ", paste(names(object@performance), collapse = ', '), ".\n", sep = '')
            else
              cat("Performance Measures: None calculated yet.\n", sep = '')
          })

setGeneric("predictions", function(object, ...)
{standardGeneric("predictions")})
setMethod("predictions", c("ClassifyResult"),
          function(object)
          {
            object@predictions
          })

setGeneric("features", function(object, ...)
{standardGeneric("features")})
setMethod("features", c("ClassifyResult"),
          function(object)
          {
            object@chosenFeatures
          })

setGeneric("performance", function(object, ...)
{standardGeneric("performance")})
setMethod("performance", c("ClassifyResult"),
          function(object)
          {
            object@performance
          })

setMethod("sampleNames", c("ClassifyResult"),
          function(object)
          {
            object@originalNames
          })

setMethod("featureNames", c("ClassifyResult"),
          function(object)
          {
            object@originalFeatures
          })

setGeneric("actualClasses", function(object, ...)
{standardGeneric("actualClasses")})
setMethod("actualClasses", c("ClassifyResult"),
          function(object)
          {
            object@actualClasses
          })