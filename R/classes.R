setClass("TransformParams", representation(
  transform = "function",
  otherParams = "list")
)

setGeneric("TransformParams", function(transform, ...)
{standardGeneric("TransformParams")})

setMethod("TransformParams", c("function"),
          function(transform, ...)
          {
            new("TransformParams", transform = transform, otherParams = list(...))
          })

setClassUnion("functionOrList", c("function", "list"))

setClass("SelectionParams", representation(
  featureSelection = "functionOrList",
  minPresence = "numeric",
  otherParams = "list")
)

setGeneric("SelectionParams", function(featureSelection, ...)
{standardGeneric("SelectionParams")})
setMethod("SelectionParams", character(0), function()
{
  new("SelectionParams", featureSelection = limmaSelection, minPresence = 1,
      otherParams = list(nFeatures = seq(100, 500, 100)))
})
setMethod("SelectionParams", c("functionOrList"),
          function(featureSelection, minPresence = 1, ...)
          {
            new("SelectionParams", featureSelection = featureSelection,
                minPresence = minPresence, otherParams = list(...))
          })

setClass("TrainParams", representation(
  classifier = "function",
  transposeExpression = "logical",
  doesTests = "logical",
  otherParams = "list")
)

setGeneric("TrainParams", function(classifier, transposeExpression, doesTests, ...)
{standardGeneric("TrainParams")})
setMethod("TrainParams", character(0), function()
{
  new("TrainParams", classifier = dlda, transposeExpression = TRUE, doesTests = FALSE)
})
setMethod("TrainParams", c("function"),
          function(classifier, transposeExpression, doesTests, ...)
          {
            new("TrainParams", classifier = classifier, transposeExpression = transposeExpression,
                doesTests = doesTests, otherParams = list(...))
          })

setClass("PredictParams", representation(
  predictor = "function",
  transposeExpression = "logical",
  multipleResults = "logical",
  getClasses = "function",
  otherParams = "list")
)

setGeneric("PredictParams", function(predictor, transposeExpression, multipleResults, getClasses, ...)
{standardGeneric("PredictParams")})
setMethod("PredictParams", character(0), function()
{
  new("PredictParams", predictor = predict, transposeExpression = TRUE, multipleResults = FALSE,
      getClasses = function(result){result[["class"]]})
})
setMethod("PredictParams", c("function"),
          function(predictor, transposeExpression, multipleResults, getClasses, ...)
          {
            new("PredictParams", predictor = predictor, transposeExpression = transposeExpression,
                multipleResults = multipleResults, getClasses = getClasses, otherParams = list(...))
          })

setGeneric("ClassifyResult", function(originalNames, originalFeatures, chosenFeatures, ...)
{standardGeneric("ClassifyResult")})
setClass("ClassifyResult", representation(
  originalNames = "character",
  originalFeatures = "character",
  chosenFeatures = "list",  
  actualClasses = "factor",
  predictions = "list",
  validation = "list",  
  errors = "list")
)
setMethod("ClassifyResult", c("character", "character", "list"),
          function(originalNames, originalFeatures, chosenFeatures, predictions, actualClasses, validation)
          {
            new("ClassifyResult", predictions = predictions, chosenFeatures = chosenFeatures,
                actualClasses = actualClasses, validation = validation,
                originalNames = originalNames, originalFeatures = originalFeatures)
          })
setMethod("show", c("ClassifyResult"),
          function(object)
          {
            cat("An object of class 'ClassifyResult'.\n")
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
            if(length(object@errors) > 0)
              cat("Errors: ", paste(names(object@errors), collapse = ', '), ".\n", sep = '')
            else
              cat("Errors: None calculated yet.\n", sep = '')
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

setGeneric("errors", function(object, ...)
{standardGeneric("errors")})
setMethod("errors", c("ClassifyResult"),
          function(object)
          {
            object@errors
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