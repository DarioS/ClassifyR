# Delete when sparsediscrim is restored to CRAN.
dlda <- function(x, ...) {
  UseMethod("dlda")
}

setOldClass("pamrtrained")
setOldClass("dlda")
setOldClass("svm")
setOldClass("multnet")
setOldClass("randomForest")

setClassUnion("functionOrNULL", c("function", "NULL"))
setClassUnion("functionOrList", c("function", "list"))
setClassUnion("integerOrNumeric", c("integer", "numeric"))
setClassUnion("characterOrDataFrame", c("character", "DataFrame"))

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

setClass("FeatureSetCollection", representation(
  sets = "list")
)

setGeneric("FeatureSetCollection", function(sets, ...)
{standardGeneric("FeatureSetCollection")})
setMethod("FeatureSetCollection", c("list"),
          function(sets)
          {
            new("FeatureSetCollection", sets = sets)
          })

setMethod("length", c("FeatureSetCollection"),
    function(x)
{
    length(x@sets)
})

setMethod("show", c("FeatureSetCollection"),
          function(object)
          {
            setType <- ifelse(class(object@sets[[1]]) == "character", "features", "interactors")
            setTypeText <- ifelse(setType == "features", "feature sets.\n", "sets of binary interactions.\n")
            if(setType == "features")
              setlElementsFunction <- length
            else
              setlElementsFunction <- nrow
            cat("An object of class 'FeatureSetCollection' consisting of", length(object@sets), setTypeText)
            setSizes <- sapply(object@sets, setlElementsFunction)
            featureText <- ifelse(setType == "features", "features.", "binary interactions.")
            cat("Smallest set:", min(setSizes), featureText, "Largest set:", max(setSizes), featureText, "\n")
            
            maxIndex <- min(length(object@sets), 3)
            featuresConcatenated <- sapply(object@sets[1:maxIndex], function(set)
            {
              if(setType == "features")
              {
                if(length(set) > 5)
                {
                  set <- set[1:6]
                  set[6] <- "..."
                }
                setText <- paste(set, collapse = ", ")
              } else { # Binary interactors
                setText <- paste(set[, 1], set[, 2], sep = '-', collapse = ", ")
                if(nrow(set) > 5)
                {
                  setText <- paste(setText, ", ...", sep = '')
                }
              }
              setText
            })
            setsText <- paste(names(object@sets)[1:maxIndex], featuresConcatenated, sep = ": ")
            setsText <- paste(setsText, collapse = '\n')
            cat(setsText)
            if(length(object@sets) > 6)
              cat("\n ...                ...\n")
            else
              cat("\n")
            minIndex <- max(length(object@sets) - 2, maxIndex + 1)
            if(minIndex <= length(object@sets))
            {
              lastIndex <- length(object@sets)
              featuresConcatenated <- sapply(object@sets[minIndex:lastIndex], function(set)
              {
                if(setType == "features")
                {
                  if(length(set) > 5)
                  {
                    set <- set[1:6]
                    set[6] <- "..."
                  }
                  setText <- paste(set, collapse = ", ")
                } else { # Binary interactors
                  setText <- paste(set[, 1], set[, 2], sep = '-', collapse = ", ")
                  if(nrow(set) > 5)
                  {
                    setText <- paste(setText, ", ...", sep = '')
                  }
                }
                setText
              })
              setsText <- paste(names(object@sets)[minIndex:lastIndex], featuresConcatenated, sep = ": ")
              setsText <- paste(setsText, collapse = '\n')
              cat(setsText)
            }
          }
)

setMethod("[", c("FeatureSetCollection", "integerOrNumeric", "missing", "ANY"),
    function(x, i, j, ..., drop = TRUE)
{
    new("FeatureSetCollection", sets = x@sets[i])
})

setMethod("[[", c("FeatureSetCollection", "ANY", "missing"),
    function(x, i, j, ...)
{
    x@sets[[i]]
})

setClassUnion("FeatureSetCollectionOrNULL", c("FeatureSetCollection", "NULL"))

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
  new("ResubstituteParams", nFeatures = seq(10, 100, 10), performanceType = "balanced error",
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
  subsetToSelections = "logical",  
  otherParams = "list")
)

setGeneric("SelectParams", function(featureSelection, ...)
{standardGeneric("SelectParams")})
setMethod("SelectParams", character(0), function()
{
  new("SelectParams", featureSelection = differentMeansSelection,
      selectionName = "Difference in Means", minPresence = 1,
      intermediate = character(0), subsetToSelections = TRUE,
      otherParams = list(resubstituteParams = ResubstituteParams()))
})
setMethod("SelectParams", c("functionOrList"),
          function(featureSelection, selectionName, minPresence = 1, intermediate = character(0),
                   subsetToSelections = TRUE, ...)
          {
            if(missing(selectionName) && !is.list(featureSelection))
              selectionName <- .methodFormals(featureSelection)[["selectionName"]]
            others <- list(...)
            if(is.list(featureSelection))
              others <- unlist(others, recursive = FALSE)
            if(is.null(others)) others <- list()
            new("SelectParams", featureSelection = featureSelection,
                selectionName = selectionName, minPresence = minPresence,
                intermediate = intermediate, subsetToSelections = subsetToSelections,
                otherParams = others)
          })

setClass("TrainParams", representation(
  classifier = "function",
  intermediate = "character",
  getFeatures = "functionOrNULL",
  otherParams = "list")
)

setGeneric("TrainParams", function(classifier, ...)
{standardGeneric("TrainParams")})
setMethod("TrainParams", character(0), function()
{
  new("TrainParams", classifier = DLDAtrainInterface, intermediate = character(0), getFeatures = NULL)
})
setMethod("TrainParams", c("function"),
          function(classifier, intermediate = character(0), getFeatures = NULL, ...)
          {
            new("TrainParams", classifier = classifier, intermediate = intermediate,
                getFeatures = getFeatures, otherParams = list(...))
          })

setClass("PredictParams", representation(
  predictor = "functionOrNULL",
  intermediate = "character",
  otherParams = "list")
)

setGeneric("PredictParams", function(predictor, ...)
{standardGeneric("PredictParams")})
setMethod("PredictParams", character(0), function()
{
  new("PredictParams", predictor = DLDApredictInterface, intermediate = character(0))
})
setMethod("PredictParams", c("functionOrNULL"),
          function(predictor, intermediate = character(0), ...)
          {
            if(missing(predictor))
              stop("Either a function or NULL must be specified by 'predictor'.")
            
            new("PredictParams", predictor = predictor, intermediate = intermediate,
                otherParams = list(...))
          })

setGeneric("SelectResult", function(dataset, selection, totalFeatures, rankedFeatures, chosenFeatures, ...)
{standardGeneric("SelectResult")})
setClass("SelectResult", representation(
  datasetName = "character",
  selectionName = "character",
  totalFeatures = "numeric",
  rankedFeatures = "list",
  chosenFeatures = "list"
))

setMethod("SelectResult", c("character", "character", "numeric", "list", "list"),
          function(dataset, selection, totalFeatures, rankedFeatures, chosenFeatures)
          {
            new("SelectResult", datasetName = dataset, selectionName = selection, 
                totalFeatures = totalFeatures, rankedFeatures = rankedFeatures,
                chosenFeatures = chosenFeatures)
          })

setMethod("show", c("SelectResult"),
          function(object)
          {
            if(class(object@chosenFeatures[[1]]) == "list")
            {
              selectionSizes <- unlist(lapply(object@chosenFeatures, function(resampling)
                                       {
                                         lapply(resampling, function(fold)
                                         {
                                               if(is.vector(fold)) length(fold)
                                               else nrow(fold) # Stored in a data frame.
                                         })
                                       })
                                       )
            } else {
              if(is.vector(object@chosenFeatures[[1]]) || "Pairs" %in% class(object@chosenFeatures[[1]]))
                selectionSizes <- sapply(object@chosenFeatures, length)
              else
                selectionSizes <- sapply(object@chosenFeatures, nrow) # Stored in a data frame.
            }
            cat("An object of class 'SelectResult'.\n")
            cat("Data Set Name: ", object@datasetName, ".\n", sep = '')
            cat("Feature Selection Name: ", object@selectionName, ".\n", sep = '')
            cat("Features Considered: ", object@totalFeatures, ".\n", sep = '')
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
  originalFeatures = "characterOrDataFrame",
  selectResult = "SelectResult",
  actualClasses = "factor",
  models = "list",
  predictions = "list",
  validation = "list",  
  performance = "list",
  tune = "list")
)
setMethod("ClassifyResult", c("character", "character", "character", "character", "character"),
          function(datasetName, classificationName, selectionName, originalNames, originalFeatures, totalFeatures,
                   rankedFeatures, chosenFeatures, models, predictions, actualClasses, validation, tune = list(NULL))
          {
            new("ClassifyResult", datasetName = datasetName, classificationName = classificationName,
                predictions = predictions, selectResult = SelectResult(datasetName, selectionName, totalFeatures, rankedFeatures, chosenFeatures),
                actualClasses = actualClasses, models = models, validation = validation,
                originalNames = originalNames, originalFeatures = originalFeatures, tune = tune)
          })
setMethod("show", c("ClassifyResult"),
          function(object)
          {
            cat("An object of class 'ClassifyResult'.\n")
            cat("Data Set Name: ", object@datasetName, ".\n", sep = '')
            cat("Classification Name: ", object@classificationName, ".\n", sep = '')
            cat("Feature Selection Name: ", object@selectResult@selectionName, ".\n", sep = '')
            if(length(unlist(object@selectResult@chosenFeatures)) == 0)
            {
              cat("Features: All used.\n")
            } else if(object@validation[[1]] != "permuteFold")
            {
              cat("Features: List of length ", length(object@selectResult@chosenFeatures), " of feature identifiers.\n", sep = '')
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
                  subListText, " of feature identifiers.\n", sep = '')
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

setGeneric("sampleNames", function(object, ...)
{standardGeneric("sampleNames")})
setMethod("sampleNames", c("ClassifyResult"),
          function(object)
          {
            object@originalNames
          })

setGeneric("featureNames", function(object, ...)
{standardGeneric("featureNames")})
setMethod("featureNames", c("ClassifyResult"),
          function(object)
          {
            object@originalFeatures
          })

setGeneric("features", function(object, ...)
{standardGeneric("features")})
setMethod("features", c("ClassifyResult"),
          function(object)
          {
            object@selectResult@chosenFeatures
          })

setGeneric("models", function(object, ...)
{standardGeneric("models")})
setMethod("models", c("ClassifyResult"),
          function(object)
          {
            object@models
          })

setGeneric("predictions", function(object, ...)
{standardGeneric("predictions")})
setMethod("predictions", c("ClassifyResult"),
          function(object)
          {
            object@predictions
          })

setGeneric("performance", function(object, ...)
{standardGeneric("performance")})
setMethod("performance", c("ClassifyResult"),
          function(object)
          {
            object@performance
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