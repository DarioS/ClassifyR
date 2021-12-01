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
setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClassUnion("characterOrMissing", c("character", "missing"))
setClassUnion("numericOrMissing", c("numeric", "missing"))
setClassUnion("BiocParallelParamOrMissing", c("BiocParallelParam", "missing"))
setClassUnion("characterOrDataFrame", c("character", "DataFrame"))
setClassUnion("listOrNULL", c("list", "NULL"))

setClass("CrossValParams", representation(
  samplesSplits = "character",
  permutations = "numericOrNULL",
  percentTest = "numericOrNULL",  
  folds = "numericOrNULL",
  leave = "numericOrNULL",
  tuneMode = "character",
  parallelParams = "BiocParallelParam",
  seed = "numeric"
  )
)

setGeneric("CrossValParams", function(samplesSplits = c("Permute k-Fold", "Permute Percentage Split", "Leave-k-Out", "k-Fold"),
                                      permutations, percentTest, folds, leave, tuneMode = c("Resubstitution", "Nested CV"),
                                      parallelParams, seed)
  standardGeneric("CrossValParams"))

setMethod("CrossValParams", c("missing", "missing", "missing", "missing",
                              "missing", "missing", "missing", "missing"),
          function()
{
  new("CrossValParams", samplesSplits = "Permute k-Fold", permutations = 100,
      folds = 5L, tuneMode = "Resubstitution", parallelParams = bpparam(),
      seed = 987654321)
})

setMethod("CrossValParams", c("characterOrMissing", "numericOrMissing", "numericOrMissing",
                              "numericOrMissing", "numericOrMissing", "characterOrMissing", "BiocParallelParamOrMissing", "numericOrMissing"),
          function(samplesSplits = c("Permute k-Fold", "Permute Percentage Split", "Leave-k-Out", "k-Fold"),
                   permutations = 100, percentTest = 25, folds = 5, leave = 2,
                   tuneMode = c("Resubstitution", "Nested CV"), parallelParams = bpparam(), seed = 987654321)
          {
            samplesSplits <- match.arg(samplesSplits)
            tuneMode <- match.arg(tuneMode)
            if(samplesSplits %in% c("Permute k-Fold", "Permute Percentage Split"))
            {
              if(missing(permutations))
                permutations <- 100
              if(samplesSplits == "Permute k-Fold")
              {
                percentTest <- NULL
                if(missing(folds))
                  folds <- 5
              } else if(samplesSplits == "Permute Percentage Split"){
                folds <- NULL
                if(missing(percentTest))
                  percentTest <- 25
              }
              leave <- NULL
            } else if(samplesSplits == "k-Fold") {
              if(missing(folds))
                folds <- 5
              permutations <- NULL
              leave <- NULL
              percentTest <- NULL
            } else if(samplesSplits == "Leave-k-Out") {
              if(missing(leave))
                 leave <- 2L
              permutations <- NULL
              percentTest <- NULL
              folds <- NULL
            }
            if(missing(parallelParams)) parallelParams <- bpparam()
            if(missing(seed)) seed <- 987654321
            new("CrossValParams", samplesSplits = samplesSplits, permutations = permutations,
                percentTest = percentTest, folds = folds, leave = leave, tuneMode = tuneMode,
                parallelParams = parallelParams, seed = seed)
          })

setClass("StageParams", representation("VIRTUAL"))
setClassUnion("StageParamsOrMissing", c("StageParams", "missing"))
setClassUnion("StageParamsOrMissingOrNULL", c("StageParams", "missing", "NULL"))

setClass("TransformParams", representation(
  transform = "function",
  characteristics = "DataFrame",
  intermediate = "character",  
  otherParams = "list"), contains = "StageParams"
)

setClassUnion("TransformParamsOrNULL", c("TransformParams", "NULL"))

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

setMethod("[", c("FeatureSetCollection", "numeric", "missing", "ANY"),
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

setClass("SelectParams", representation(
  featureRanking = "functionOrList",
  characteristics = "DataFrame",
  minPresence = "numeric",
  intermediate = "character",
  subsetToSelections = "logical",
  tuneParams = "listOrNULL",
  otherParams = "listOrNULL"), contains = "StageParams"
)

setClassUnion("SelectParamsOrNULL", c("SelectParams", "NULL"))

setGeneric("SelectParams", function(featureRanking, ...)
standardGeneric("SelectParams"))
setMethod("SelectParams", "missing", function()
{
  new("SelectParams", featureRanking = differentMeansRanking,
      characteristics = DataFrame(characteristic = "Selection Name", value = "Difference in Means"),
      minPresence = 1, intermediate = character(0), subsetToSelections = TRUE,
      tuneParams = list(nFeatures = seq(10, 100, 10), performanceType = "Balanced Error"))
})
setMethod("SelectParams", c("functionOrList"),
          function(featureRanking, characteristics = DataFrame(), minPresence = 1, 
                   intermediate = character(0), subsetToSelections = TRUE, tuneParams = list(nFeatures = seq(10, 100, 10), performanceType = "Balanced Error"), ...)
          {
            if(!is.list(featureRanking) && (ncol(characteristics) == 0 || !"Selection Name" %in% characteristics[, "characteristic"]))
            {
              characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Selection Name", value = .ClassifyRenvir[["functionsTable"]][.ClassifyRenvir[["functionsTable"]][, "character"] == featureRanking@generic, "name"]))
            }
            if(is.list(featureRanking) && (ncol(characteristics) == 0 || !"Ensemble Selection" %in% characteristics[, "characteristic"]))
            {
              selectMethodNames <- unlist(lapply(featureRanking, function(rankingFunction) .ClassifyRenvir[["functionsTable"]][.ClassifyRenvir[["functionsTable"]][, "character"] == rankingFunction@generic, "name"]))
              characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Ensemble Selection", value = paste(selectMethodNames, collapse = ", ")))
            }
            others <- list(...)
            if(length(others) == 0) others <- NULL
            new("SelectParams", featureRanking = featureRanking,
                characteristics = characteristics, minPresence = minPresence,
                intermediate = intermediate, subsetToSelections = subsetToSelections,
                tuneParams = tuneParams, otherParams = others)
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
  tuneParams = "listOrNULL",
  otherParams = "listOrNULL",
  getFeatures = "functionOrNULL"
), contains = "StageParams"
)

setGeneric("TrainParams", function(classifier, ...) standardGeneric("TrainParams"))
setMethod("TrainParams", "missing", function()
{
  new("TrainParams", classifier = DLDAtrainInterface,
      characteristics = S4Vectors::DataFrame(characteristic = "Classifier Name", value = "Diagonal LDA"),
      intermediate = character(0), getFeatures = NULL)
})
setMethod("TrainParams", c("function"),
          function(classifier, balancing = c("downsample", "upsample", "none"), characteristics = DataFrame(), intermediate = character(0), tuneParams = NULL, getFeatures = NULL, ...)
          {
            if(ncol(characteristics) == 0 || !"Classifier Name" %in% characteristics[, "characteristic"])
            {
              characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Classifier Name", value = .ClassifyRenvir[["functionsTable"]][.ClassifyRenvir[["functionsTable"]][, "character"] == classifier@generic, "name"]))
            }
            new("TrainParams", classifier = classifier, characteristics = characteristics,
                intermediate = intermediate, getFeatures = getFeatures, tuneParams = tuneParams,
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
  otherParams = "listOrNULL"), contains = "StageParams"
)

setGeneric("PredictParams", function(predictor, ...)
standardGeneric("PredictParams"))
setMethod("PredictParams", "missing", function()
{
  new("PredictParams", predictor = DLDApredictInterface,
      characteristics = S4Vectors::DataFrame(characteristic = "Predictor Name", value = "Diagonal LDA"),
      intermediate = character(0), otherParams = NULL)
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
            others <- list(...)
            if(length(others) == 0) others <- NULL
            new("PredictParams", predictor = predictor, characteristics = characteristics,
                intermediate = intermediate, otherParams = others)
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

setClassUnion("PredictParamsOrNULL", c("PredictParams", "NULL"))

setGeneric("ModellingParams", function(balancing, transformParams, selectParams, trainParams, predictParams)
  standardGeneric("ModellingParams"))
setClass("ModellingParams", representation(
  balancing = "character",
  transformParams = "TransformParamsOrNULL",
  selectParams = "SelectParamsOrNULL",
  trainParams = "TrainParams",
  predictParams = "PredictParamsOrNULL"
))
setMethod("ModellingParams", c("missing", "missing", "missing", "missing", "missing"), function()
{
  new("ModellingParams", balancing = "downsample", transformParams = NULL, selectParams = SelectParams(),
      trainParams = TrainParams(), predictParams = PredictParams())
})
setMethod("ModellingParams", c("characterOrMissing", "StageParamsOrMissingOrNULL", "StageParamsOrMissingOrNULL", "StageParamsOrMissing", "StageParamsOrMissingOrNULL"),
          function(balancing = c("downsample", "upsample", "none"),
                   transformParams = NULL, selectParams = SelectParams(),
                   trainParams = TrainParams(), predictParams = PredictParams())
          {
            if(missing(balancing)) balancing <- "downsample"
            balancing <- match.arg(balancing)
            if(missing(transformParams)) transformParams <- NULL
            if(missing(selectParams)) selectParams <- SelectParams()
            if(missing(trainParams)) trainParams <- TrainParams()
            if(missing(predictParams)) predictParams <- PredictParams()
            
            new("ModellingParams", balancing = balancing, transformParams = transformParams,
                selectParams = selectParams, trainParams = trainParams, predictParams = predictParams)
          })

setGeneric("ClassifyResult", function(characteristics, originalNames, originalFeatures, ...)
standardGeneric("ClassifyResult"))
setClass("ClassifyResult", representation(
  characteristics = "DataFrame",
  originalNames = "character",
  originalFeatures = "characterOrDataFrame",
  rankedFeatures = "listOrNULL",
  chosenFeatures = "listOrNULL",
  actualClasses = "factor",
  models = "list",
  tune = "listOrNULL",
  predictions = "data.frame",
  performance = "listOrNULL")
)
setMethod("ClassifyResult", c("DataFrame", "character", "characterOrDataFrame"),
          function(characteristics, originalNames, originalFeatures,
                   rankedFeatures, chosenFeatures, models, tunedParameters, predictions, actualClasses)
          {
            new("ClassifyResult", characteristics = characteristics,
                originalNames = originalNames, originalFeatures = originalFeatures,
                rankedFeatures = rankedFeatures, chosenFeatures = chosenFeatures,
                models = models, tune = tunedParameters,
                predictions = predictions, actualClasses = actualClasses)
          })
setMethod("show", "ClassifyResult", function(object)
          {
            cat("An object of class 'ClassifyResult'.\n")
            cat("Characteristics:\n")
            print(as.data.frame(object@characteristics), row.names = FALSE)
            
            cat("Features: List of length ", length(object@chosenFeatures), " of feature identifiers.\n", sep = '')
            cat("Predictions: A data frame of ", nrow(object@predictions), " rows.\n", sep = '')
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
            object@chosenFeatures
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
              nrow(predictions(result))
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

