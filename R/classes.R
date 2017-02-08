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

setClass("SelectParams", representation(
  featureSelection = "functionOrList",
  selectionName = "character",
  minPresence = "numeric",
  intermediate = "character",
  subsetExpressionData = "logical",  
  otherParams = "list")
)

setGeneric("SelectParams", function(featureSelection, ...)
{standardGeneric("SelectParams")})
setMethod("SelectParams", character(0), function()
{
  new("SelectParams", featureSelection = limmaSelection,
      selectionName = "Limma moderated t-test", minPresence = 1,
      intermediate = character(0), subsetExpressionData = TRUE,
      otherParams = list(resubstituteParams = ResubstituteParams()))
})
setMethod("SelectParams", c("functionOrList"),
          function(featureSelection, selectionName, minPresence = 1, intermediate = character(0),
                   subsetExpressionData = TRUE, ...)
          {
            if(missing(selectionName) && !is.list(featureSelection))
              selectionName <- .methodFormals(featureSelection, "ExpressionSet")[["selectionName"]]
            others <- list(...)
            if(is.list(featureSelection))
              others <- unlist(others, recursive = FALSE)
            if(is.null(others)) others <- list()
            new("SelectParams", featureSelection = featureSelection,
                selectionName = selectionName, minPresence = minPresence,
                intermediate = intermediate, subsetExpressionData = subsetExpressionData,
                otherParams = others)
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
            new("TrainParams", classifier = classifier,
                transposeExpression = transposeExpression, doesTests = doesTests,
                intermediate = intermediate, otherParams = list(...))
          })

setClass("PredictParams", representation(
  predictor = "function",
  transposeExpression = "logical",
  intermediate = "character",    
  getClasses = "function",
  otherParams = "list")
)

setGeneric("PredictParams", function(predictor, ...)
{standardGeneric("PredictParams")})
setMethod("PredictParams", character(0), function()
{
  new("PredictParams", predictor = predict, transposeExpression = TRUE,
      intermediate = character(0), getClasses = function(result){result[["class"]]})
})
setMethod("PredictParams", c("function"),
          function(predictor, transposeExpression, intermediate = character(0), getClasses, ...)
          {
            new("PredictParams", predictor = predictor, transposeExpression = transposeExpression,
                intermediate = intermediate, getClasses = getClasses, otherParams = list(...))
          })

setGeneric("SelectResult", function(dataset, selection, rankedFeatures, chosenFeatures, ...)
{standardGeneric("SelectResult")})
setClass("SelectResult", representation(
  datasetName = "character",
  selectionName = "character",
  rankedFeatures = "list",
  chosenFeatures = "list"
))

setMethod("SelectResult", c("character", "character", "list", "list"),
          function(dataset, selection, rankedFeatures, chosenFeatures)
          {
            new("SelectResult", datasetName = dataset, selectionName = selection, 
                rankedFeatures = rankedFeatures, chosenFeatures = chosenFeatures)
          })
setMethod("show", c("SelectResult"),
          function(object)
          {
            if(class(object@chosenFeatures[[1]]) == "list")
            {
              selectionSizes <- unlist(lapply(object@chosenFeatures, function(resampling)
                lapply(resampling, function(fold) length(fold))))
            } else {selectionSizes <- sapply(object@chosenFeatures, length)}
            cat("An object of class 'SelectResult'.\n")
            cat("Dataset Name: ", object@datasetName, ".\n", sep = '')
            cat("Feature Selection Name: ", object@selectionName, ".\n", sep = '')
            if(length(object@rankedFeatures) > 0) # Some methods don't use rankings.
            {
              if(class(object@rankedFeatures[[1]]) == "list") # Must be from resampling and folding.
                featureLength <- length(object@rankedFeatures[[1]][[1]])
              else # Either from direct feature selection, or a cross-validation method that doesn't have folds.
                featureLength <- length(object@rankedFeatures[[1]])
              cat("Features Considered: ", featureLength, ".\n", sep = '')
            }
            selectionsText <- paste("Selections: List of length", length(object@chosenFeatures))
            if(class(object@chosenFeatures[[1]]) == "list")
              selectionsText <- paste(selectionsText, "of lists of length", length(object@chosenFeatures[[1]]))
            cat(selectionsText, ".\n", sep = '')
            if(length(selectionSizes) > 1)
              cat("Selection Size Range: Between ", min(selectionSizes), " and ", max(selectionSizes), " features.\n", sep = '')
            else
              cat("Selection Size: ", selectionSizes[[1]], " features.\n", sep = '')
          })

setGeneric("ClassifyResult", function(datasetName, classificationName, selectionName, originalNames, originalFeatures, ...)
{standardGeneric("ClassifyResult")})
setClass("ClassifyResult", representation(
  datasetName = "character",
  classificationName = "character",
  originalNames = "character",
  originalFeatures = "character",
  selectResult = "SelectResult",
  actualClasses = "factor",
  predictions = "list",
  validation = "list",  
  performance = "list",
  tune = "list")
)
setMethod("ClassifyResult", c("character", "character", "character", "character", "character"),
          function(datasetName, classificationName, selectionName, originalNames, originalFeatures,
                   rankedFeatures, chosenFeatures, predictions, actualClasses, validation, tune = list(NULL))
          {
            new("ClassifyResult", datasetName = datasetName, classificationName = classificationName,
                predictions = predictions, selectResult = SelectResult(datasetName, selectionName, rankedFeatures, chosenFeatures),
                actualClasses = actualClasses, validation = validation,
                originalNames = originalNames, originalFeatures = originalFeatures, tune = tune)
          })
setMethod("show", c("ClassifyResult"),
          function(object)
          {
            cat("An object of class 'ClassifyResult'.\n")
            cat("Dataset Name: ", object@datasetName, ".\n", sep = '')
            cat("Classification Name: ", object@classificationName, ".\n", sep = '')
            cat("Feature Selection Name: ", object@selectResult@selectionName, ".\n", sep = '')
            if(object@validation[[1]] != "resampleFold")
            {
              cat("Features: List of length ", length(object@selectResult@chosenFeatures), " of row indices.\n", sep = '')
            } else # Resample and fold. Nested lists.
            {
              elementsLengths <- sapply(object@selectResult@chosenFeatures, length)
              if(diff(range(elementsLengths)) == 0)
              {
                subListText <- paste("length", unique(elementsLengths))
              } else
              {
                subListText <- paste("lengths between", min(elementsLengths), "and", max(elementsLengths))
              }
              cat("Features: List of length ", length(object@selectResult@chosenFeatures), " of lists of ",
                  subListText, " of row indices.\n", sep = '')
            }            
            cat("Validation: ")
            cat(.validationText(object), ".\n", sep = '')
            cat("Predictions: List of data frames of length ", length(object@predictions),
                ".\n", sep = '')
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
            object@selectResult@chosenFeatures
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

setGeneric("tunedParameters", function(object, ...)
{standardGeneric("tunedParameters")})
setMethod("tunedParameters", c("ClassifyResult"),
          function(object)
          {
            object@tune
          })

setGeneric("totalPredictions", function(result, ...)
{standardGeneric("totalPredictions")})
setMethod("totalPredictions", c("ClassifyResult"),
          function(result)
          {
              nrow(do.call(rbind, predictions(result)))
          })