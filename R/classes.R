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
setClassUnion("listOrNULL", c("list", "NULL"))
setClassUnion("listOrCharacterOrNULL", c("list", "character", "NULL"))

setClass("TransformParams", representation(
  transform = "function",
  characteristics = "DataFrame",
  intermediate = "character",  
  otherParams = "list")
)

setGeneric("TransformParams", function(transform, ...)
standardGeneric("TransformParams"))

setMethod("TransformParams", c("function"),
          function(transform, characteristics = DataFrame(), intermediate = character(0), ...)
          {
            if(ncol(characteristics) == 0 || !"Transform Name" %in% characteristics[, "characteristic"])
            {
              characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Transform Name", value = .ClassifyRenvir[["functionsTable"]][.ClassifyRenvir[["functionsTable"]][, "character"] == transform@generic, "name"]))
            }
            new("TransformParams", transform = transform, characteristics = characteristics,
                intermediate = intermediate, otherParams = list(...))
          })

setMethod("show", "TransformParams",
          function(object)
          {
            cat("An object of class 'TransformParams'.\n")
            index <- na.omit(match("Transform Name", object@characteristics[, "characteristic"]))
            transText <- object@characteristics[index, "characteristic"]
            cat("Transform Name: ", object@characteristics[index, "value"], ".\n", sep = '')
            
            otherInfo <- object@characteristics[-index, ]
            if(nrow(otherInfo) > 0)
            {
              for(rowIndex in 1:nrow(otherInfo))
              {
                cat(otherInfo[rowIndex, "characteristic"], ": ", otherInfo[rowIndex, "value"], ".\n", sep = '')
              }
            }
          })

setClass("FeatureSetCollection", representation(
  sets = "list")
)

setGeneric("FeatureSetCollection", function(sets, ...)
standardGeneric("FeatureSetCollection"))
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

setMethod("show", "FeatureSetCollection",
          function(object)
          {
            setType <- ifelse("character" %in% class(object@sets[[1]]), "features", "interactors")
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
  otherParams = "list")
)

setGeneric("ResubstituteParams", function(nFeatures, performanceType, ...)
standardGeneric("ResubstituteParams"))

setMethod("ResubstituteParams", "missing", function()
{
  new("ResubstituteParams", nFeatures = seq(10, 100, 10), performanceType = "balanced error")
})

setMethod("ResubstituteParams", c("numeric", "character"),
          function(nFeatures, performanceType, ...)
          {
            new("ResubstituteParams", nFeatures = nFeatures, performanceType = performanceType,
                otherParams = list(...))
          })

setClass("SelectParams", representation(
  featureSelection = "functionOrList",
  characteristics = "DataFrame",
  minPresence = "numeric",
  intermediate = "character",
  subsetToSelections = "logical",  
  otherParams = "list")
)

setGeneric("SelectParams", function(featureSelection, ...)
standardGeneric("SelectParams"))
setMethod("SelectParams", "missing", function()
{
  new("SelectParams", featureSelection = differentMeansSelection,
      characteristics = DataFrame(characteristic = "Selection Name", value = "Difference in Means"),
      minPresence = 1, intermediate = character(0), subsetToSelections = TRUE,
      otherParams = list(resubstituteParams = ResubstituteParams()))
})
setMethod("SelectParams", c("functionOrList"),
          function(featureSelection, characteristics = DataFrame(), minPresence = 1, 
                   intermediate = character(0), subsetToSelections = TRUE, ...)
          {
            if(!is.list(featureSelection) && (ncol(characteristics) == 0 || !"Selection Name" %in% characteristics[, "characteristic"]))
            {
              characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Selection Name", value = .ClassifyRenvir[["functionsTable"]][.ClassifyRenvir[["functionsTable"]][, "character"] == featureSelection@generic, "name"]))
            }
            if(is.list(featureSelection) && (ncol(characteristics) == 0 || !"Ensemble Selection" %in% characteristics[, "characteristic"]))
            {
              selectMethodNames <- unlist(lapply(featureSelection, function(selectFunction) .ClassifyRenvir[["functionsTable"]][.ClassifyRenvir[["functionsTable"]][, "character"] == selectFunction@generic, "name"]))
              characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Ensemble Selection", value = paste(selectMethodNames, collapse = ", ")))
            }
            others <- list(...)
            if(is.list(featureSelection))
              others <- unlist(others, recursive = FALSE)
            if(is.null(others)) others <- list()
            new("SelectParams", featureSelection = featureSelection,
                characteristics = characteristics, minPresence = minPresence,
                intermediate = intermediate, subsetToSelections = subsetToSelections,
                otherParams = others)
          })

setMethod("show", "SelectParams",
          function(object)
          {
            cat("An object of class 'SelectParams'.\n")
            IDs <- c("Selection Name", "Ensemble Selection")
            index <- na.omit(match(IDs, object@characteristics[, "characteristic"]))
            selectText <- object@characteristics[index, "characteristic"]
            cat(selectText, ": ", object@characteristics[index, "value"], ".\n", sep = '')
            if(selectText == "Ensemble Selection")
              cat("Minimum Functions Selected By:", object@minPresence)
            
            otherInfo <- object@characteristics[-index, ]
            if(nrow(otherInfo) > 0)
            {
              for(rowIndex in 1:nrow(otherInfo))
              {
                cat(otherInfo[rowIndex, "characteristic"], ": ", otherInfo[rowIndex, "value"], ".\n", sep = '')
              }
            }
          })

setClass("TrainParams", representation(
  classifier = "function",
  characteristics = "DataFrame",
  intermediate = "character",
  getFeatures = "functionOrNULL",
  otherParams = "list")
)

setGeneric("TrainParams", function(classifier, ...) standardGeneric("TrainParams"))
setMethod("TrainParams", "missing", function()
{
  new("TrainParams", classifier = DLDAtrainInterface,
      characteristics = S4Vectors::DataFrame(characteristic = "Classifier Name", value = "Diagonal LDA"),
      intermediate = character(0), getFeatures = NULL)
})
setMethod("TrainParams", c("function"),
          function(classifier, characteristics = DataFrame(), intermediate = character(0), getFeatures = NULL, ...)
          {
            if(ncol(characteristics) == 0 || !"Classifier Name" %in% characteristics[, "characteristic"])
            {
              characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Classifier Name", value = .ClassifyRenvir[["functionsTable"]][.ClassifyRenvir[["functionsTable"]][, "character"] == classifier@generic, "name"]))
            }
            new("TrainParams", classifier = classifier, characteristics = characteristics,
                intermediate = intermediate, getFeatures = getFeatures,
                otherParams = list(...))
          })

setMethod("show", "TrainParams",
          function(object)
          {
            cat("An object of class 'TrainParams'.\n")
            index <- na.omit(match("Classifier Name", object@characteristics[, "characteristic"]))
            trainText <- object@characteristics[index, "characteristic"]
            cat("Classifier Name: ", object@characteristics[index, "value"], ".\n", sep = '')
            
            otherInfo <- object@characteristics[-index, ]
            if(nrow(otherInfo) > 0)
            {
              for(rowIndex in 1:nrow(otherInfo))
              {
                cat(otherInfo[rowIndex, "characteristic"], ": ", otherInfo[rowIndex, "value"], ".\n", sep = '')
              }
            }

            if(!is.null(object@getFeatures))
              cat("Selected Features Extracted By: ", object@getFeatures@generic, ".\n", sep = '')
          })

setClass("PredictParams", representation(
  predictor = "functionOrNULL",
  characteristics = "DataFrame",  
  intermediate = "character",
  otherParams = "list")
)

setGeneric("PredictParams", function(predictor, ...)
standardGeneric("PredictParams"))
setMethod("PredictParams", "missing", function()
{
  new("PredictParams", predictor = DLDApredictInterface,
      characteristics = S4Vectors::DataFrame(characteristic = "Predictor Name", value = "Diagonal LDA"),
      intermediate = character(0))
})
setMethod("PredictParams", c("functionOrNULL"),
          function(predictor, characteristics = DataFrame(), intermediate = character(0), ...)
          {
            if(missing(predictor))
              stop("Either a function or NULL must be specified by 'predictor'.")
            if(!is.null(predictor) && (ncol(characteristics) == 0 || !"Predictor Name" %in% characteristics[, "characteristic"]))
            {
              characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Predictor Name", value = .ClassifyRenvir[["functionsTable"]][.ClassifyRenvir[["functionsTable"]][, "character"] == predictor@generic, "name"]))
            }
            new("PredictParams", predictor = predictor, characteristics = characteristics,
                intermediate = intermediate, otherParams = list(...))
          })

setMethod("show", "PredictParams",
          function(object)
          {
            cat("An object of class 'PredictParams'.\n")
            if(ncol(object@characteristics) > 0)
            {
              index <- na.omit(match("Predictor Name", object@characteristics[, "characteristic"]))
              if(length(index) > 0)
                cat("Predictor Name: ", object@characteristics[index, "value"], ".\n", sep = '')
              
              otherInfo <- object@characteristics[-index, ]
              if(nrow(otherInfo) > 0)
              {
                for(rowIndex in 1:nrow(otherInfo))
                {
                  cat(otherInfo[rowIndex, "characteristic"], ": ", otherInfo[rowIndex, "value"], ".\n", sep = '')
                }
              }
            } 
            if(is.null(object@predictor))
              cat("Prediction is done by function specified to TrainParams.\n")
          })

setGeneric("SelectResult", function(totalFeatures, rankedFeatures, chosenFeatures, ...)
standardGeneric("SelectResult"))
setClass("SelectResult", representation(
  totalFeatures = "numeric",
  rankedFeatures = "list",
  chosenFeatures = "list"
))

setMethod("SelectResult", c("numeric", "list", "list"),
          function(totalFeatures, rankedFeatures, chosenFeatures)
          {
            new("SelectResult", totalFeatures = totalFeatures, rankedFeatures = rankedFeatures,
                chosenFeatures = chosenFeatures)
          })

setMethod("show", "SelectResult",
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

setGeneric("ClassifyResult", function(characteristics, originalNames, originalFeatures, ...)
standardGeneric("ClassifyResult"))
setClass("ClassifyResult", representation(
  characteristics = "DataFrame",
  originalNames = "character",
  originalFeatures = "characterOrDataFrame",
  selectResult = "SelectResult",
  actualClasses = "factor",
  models = "list",
  predictions = "list",
  validation = "list",  
  performance = "listOrNULL",
  tune = "listOrNULL")
)
setMethod("ClassifyResult", c("DataFrame", "character", "characterOrDataFrame"),
          function(characteristics, originalNames, originalFeatures, totalFeatures,
                   rankedFeatures, chosenFeatures, models, predictions, actualClasses, validation, tune = NULL)
          {
            new("ClassifyResult", characteristics = characteristics,
                predictions = predictions, selectResult = SelectResult(totalFeatures, rankedFeatures, chosenFeatures),
                actualClasses = actualClasses, models = models, validation = validation,
                originalNames = originalNames, originalFeatures = originalFeatures, tune = tune)
          })
setMethod("show", "ClassifyResult", function(object)
          {
            cat("An object of class 'ClassifyResult'.\n")
            cat("Characteristics:\n")
            print(as.data.frame(object@characteristics), row.names = FALSE)
            
            if(object@validation[[1]] != "permuteFold")
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
            cat("Predictions: List of data frames of length ", length(object@predictions),
                ".\n", sep = '')
            if(length(object@performance) > 0)
              cat("Performance Measures: ", paste(names(object@performance), collapse = ', '), ".\n", sep = '')
            else
              cat("Performance Measures: None calculated yet.\n", sep = '')
          })

setGeneric("sampleNames", function(object, ...)
standardGeneric("sampleNames"))
setMethod("sampleNames", c("ClassifyResult"),
          function(object)
          {
            object@originalNames
          })

setGeneric("featureNames", function(object, ...)
standardGeneric("featureNames"))
setMethod("featureNames", c("ClassifyResult"),
          function(object)
          {
            object@originalFeatures
          })

setGeneric("features", function(object, ...)
standardGeneric("features"))
setMethod("features", c("ClassifyResult"),
          function(object)
          {
            object@selectResult@chosenFeatures
          })

setGeneric("models", function(object, ...)
standardGeneric("models"))
setMethod("models", c("ClassifyResult"),
          function(object)
          {
            object@models
          })

setGeneric("predictions", function(object, ...)
standardGeneric("predictions"))
setMethod("predictions", c("ClassifyResult"),
          function(object)
          {
            object@predictions
          })

setGeneric("performance", function(object, ...)
standardGeneric("performance"))
setMethod("performance", c("ClassifyResult"),
          function(object)
          {
            object@performance
          })

setGeneric("actualClasses", function(object, ...)
standardGeneric("actualClasses"))
setMethod("actualClasses", c("ClassifyResult"),
          function(object)
          {
            object@actualClasses
          })

setGeneric("tunedParameters", function(object, ...)
standardGeneric("tunedParameters"))
setMethod("tunedParameters", "ClassifyResult",
          function(object)
          {
            object@tune
          })

setGeneric("totalPredictions", function(result, ...)
standardGeneric("totalPredictions"))
setMethod("totalPredictions", c("ClassifyResult"),
          function(result)
          {
              nrow(do.call(rbind, predictions(result)))
          })

setClass("EasyHardClassifier", representation(
  easyClassifier = "listOrNULL",
  hardClassifier = "listOrCharacterOrNULL",
  datasetIDs = "character"
))
setGeneric("EasyHardClassifier", function(easyClassifier, hardClassifier, datasetIDs)
standardGeneric("EasyHardClassifier"))
setMethod("EasyHardClassifier", c("listOrNULL", "listOrCharacterOrNULL", "character"),
          function(easyClassifier, hardClassifier, datasetIDs)
          {
            new("EasyHardClassifier", easyClassifier = easyClassifier, hardClassifier = hardClassifier,
                datasetIDs = datasetIDs)
          })

setMethod("show", "EasyHardClassifier",
          function(object)
          {
            cat("An object of class 'EasyHardClassifier'.\n")
            if(!is.null(object@easyClassifier)) easyText <- paste("A set of", length(object@easyClassifier), "rules trained on", object@datasetIDs["easy"], "data")
            else easyText <- "None"
            cat("Easy Classifier: ", easyText, ".\n", sep = '')
            if(!is.null(object@hardClassifier))
              hardText <- paste("An object of class '", class(object@hardClassifier[["model"]]), "' trained on ", object@datasetIDs["hard"], " data", sep = '')
            else hardText <- "None"
            cat("Hard Classifier: ", hardText, ".\n", sep = '')
          })

setClass("MixModelsListsSet", representation(
  set = "list")
)

setGeneric("MixModelsListsSet", function(set, ...)
standardGeneric("MixModelsListsSet"))
setMethod("MixModelsListsSet", c("list"),
          function(set)
          {
            new("MixModelsListsSet", set = set)
          })
