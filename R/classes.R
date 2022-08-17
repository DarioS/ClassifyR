##### Table of contents #####
# Set old classes
# Create union of classes
# Set generic accessors

################################################################################
#
# Set old classes
#
################################################################################



# Delete when sparsediscrim is restored to CRAN.
# Trained dlda Object
dlda <- function(x, ...) {
  UseMethod("dlda")
}
setOldClass("dlda")

# Trained pamr Object
setOldClass("pamrtrained")


# Trained svm Object
setOldClass("svm")

# Trained multnet Object
setOldClass("multnet")

# Trained coxnet Object
setOldClass("coxnet")

# Trained randomForest Object
setOldClass("randomForest")

# Trained coxph Object
setOldClass("coxph")

# Survival Data Container
setOldClass("Surv")

# Survival Forest Data Container
setOldClass("rfsrc")

################################################################################
#
# Create union of classes
#
################################################################################


# Union of A Function and NULL
setClassUnion("functionOrNULL", c("function", "NULL"))

# Union of Functions and List of Functions. Useful for allowing ensemble feature selection.
setClassUnion("functionOrList", c("function", "list"))


# Union of A Numeric Value and NULL
setClassUnion("numericOrNULL", c("numeric", "NULL"))

# Union of a Character and a DataFrame
setClassUnion("characterOrDataFrame", c("character", "DataFrame"))

# Union of a Surv class and a factor
setClassUnion("factorOrSurv", c("factor", "Surv"))

# Union of a List and NULL
setClassUnion("listOrNULL", c("list", "NULL"))

# Union of NULL and DataFrame Class
setClassUnion("DataFrameOrNULL", c("DataFrame", "NULL"))


################################################################################
#
# Params
#
################################################################################

##### CrossValParams #####

# Parameters for Cross-validation Specification

setClass("CrossValParams", representation(
    samplesSplits = "character",
    permutations = "numericOrNULL",
    percentTest = "numericOrNULL",  
    folds = "numericOrNULL",
    leave = "numericOrNULL",
    tuneMode = "character",
    adaptiveResamplingDelta = "numericOrNULL",
    parallelParams = "BiocParallelParam"
)
)

# CrossValParams constructor is an ordinary function and not S4 method for performance reasons.
CrossValParams <- function(samplesSplits = c("Permute k-Fold", "Permute Percentage Split", "Leave-k-Out", "k-Fold"),
                           permutations = 100, percentTest = 25, folds = 5, leave = 2,
                           tuneMode = c("Resubstitution", "Nested CV", "none"), adaptiveResamplingDelta = NULL, parallelParams = bpparam())
{
  samplesSplits <- match.arg(samplesSplits)
  tuneMode <- match.arg(tuneMode)
  if(samplesSplits == "Permute k-Fold")
  {
    percentTest <- NULL
    leave <- NULL
  } else if(samplesSplits == "Permute Percentage Split"){
    folds <- NULL
    leave <- NULL
  } else if(samplesSplits == "k-Fold") {
    permutations <- NULL
    leave <- NULL
    percentTest <- NULL
  } else if(samplesSplits == "Leave-k-Out") {
    permutations <- NULL
    percentTest <- NULL
    folds <- NULL
  }

  new("CrossValParams", samplesSplits = samplesSplits, permutations = permutations,
      percentTest = percentTest, folds = folds, leave = leave, tuneMode = tuneMode,
      adaptiveResamplingDelta = adaptiveResamplingDelta, parallelParams = parallelParams)
}

##### StageParams #####

# StageParams Virtual Class. A class for any one of TransformParams,
# SelectParams, TrainParams, PredictParams. Allows a method to dispatch on
# any of the parameter objects specifying any stage of cross-validation.

setClass("StageParams", representation("VIRTUAL"))
setClassUnion("StageParamsOrMissing", c("StageParams", "missing"))

# Union of a StageParams object and NULL for parameter that's optional.
setClassUnion("StageParamsOrMissingOrNULL", c("StageParams", "missing", "NULL"))


##### TransformParams #####
setClass("TransformParams", representation(
  transform = "function",
  characteristics = "DataFrame",
  intermediate = "character",  
  otherParams = "list"), contains = "StageParams"
)

# Union of a TransformParams pbject and NULL. 
setClassUnion("TransformParamsOrNULL", c("TransformParams", "NULL"))

# Parameters for Data Transformation within CV.

setGeneric("TransformParams", function(transform, ...)
standardGeneric("TransformParams"))

setMethod("TransformParams", "function",
          function(transform, characteristics = S4Vectors::DataFrame(), intermediate = character(0), ...)
          {
            if(ncol(characteristics) == 0 || !"Transform Name" %in% characteristics[, "characteristic"])
            {
              characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Transform Name", value = .ClassifyRenvir[["functionsTable"]][.ClassifyRenvir[["functionsTable"]][, "character"] == attr(transform, "name"), "name"]))
            }
            new("TransformParams", transform = transform, characteristics = characteristics,
                intermediate = intermediate, otherParams = list(...))
          })

#' Inspect Data Transformation Details
#'
#' @rdname TransformParams-class
#' @param object An object of class \code{TransformParams} to inspect.
#' @export
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

##### FeatureSetCollection #####
#' @docType class
#' @exportClass FeatureSetCollection
#' @rdname FeatureSetCollection
setClass("FeatureSetCollection", representation(sets = "list"))
#'
#' Container for Storing A Collection of Sets
#' 
#' This container is the required storage format for a collection of sets.
#' Typically, the elements of a set will either be a set of proteins (i.e.
#' character vector) which perform a particular biological process or a set of
#' binary interactions (i.e. Two-column matrix of feature identifiers).
#' 
#' 
#' @name FeatureSetCollection
#' @aliases FeatureSetCollection FeatureSetCollection-class
#' FeatureSetCollection,list-method length,FeatureSetCollection-method
#' show,FeatureSetCollection-method
#' [,FeatureSetCollection,numeric,missing,ANY-method
#' [[,FeatureSetCollection,ANY,missing-method
#' @docType class
#' @usage NULL
#' @section Constructor:
#' \describe{\item{}{
#' \code{FeatureSetCollection(sets)}}
#' }
#' \describe{
#' \item{\code{sets}}{A named list. The names of the list
#' describe the sets and the elements of the list specify the features which
#' comprise the sets.}
#' }
#' @section Summary:
#' \code{featureSets} is a \code{FeatureSetCollection} object.
#' 
#' \describe{
#'   \item{}{
#'     \code{show(featureSets)}: Prints a short summary of what \code{featureSets} contains.
#'    }
#'   \item{}{
#'     \code{length(featureSets)}: Prints how many sets of features there are.
#'    }  
#' }
#' @section Subsetting:
#' The \code{FeatureSetCollection} may be subsetted to a smaller set of elements or a single set
#' may be extracted as a vector.
#'  \code{featureSets} is a \code{FeatureSetCollection} object.
#'  \describe{
#'  \item{}{
#'    \code{featureSets[i:j]}:
#'      Reduces the object to a subset of the feature sets between elements \code{i} and \code{j}
#'      of the collection.
#'  }
#'  \item{}{
#'    \code{featureSets[[i]]}:
#'      Extract the feature set identified by \code{i}. \code{i} may be a numeric index
#'      or the character name of a feature set.
#'    }    
#'  }
#' @author Dario Strbenac
#' @examples
#' 
#'     ontology <- list(c("SESN1", "PRDX1", "PRDX2", "PRDX3", "PRDX4", "PRDX5", "PRDX6",
#'                        "LRRK2", "PARK7"),
#'                      c("ATP7A", "CCS", "NQO1", "PARK7", "SOD1", "SOD2", "SOD3",
#'                        "SZT2", "TNF"),
#'                      c("AARS", "AIMP2", "CARS", "GARS", "KARS", "NARS", "NARS2",
#'                        "LARS2", "NARS", "NARS2", "RGN", "UBA7"),
#'                      c("CRY1", "CRY2", "ONP1SW", "OPN4", "RGR"),
#'                      c("ESRRG", "RARA", "RARB", "RARG", "RXRA", "RXRB", "RXRG"),
#'                      c("CD36", "CD47", "F2", "SDC4"),
#'                      c("BUD31", "PARK7", "RWDD1", "TAF1")
#'                      )
#'     names(ontology) <- c("Peroxiredoxin Activity", "Superoxide Dismutase Activity",
#'                          "Ligase Activity", "Photoreceptor Activity",
#'                          "Retinoic Acid Receptor Activity",
#'                          "Thrombospondin Receptor Activity",
#'                          "Regulation of Androgen Receptor Activity")
#'                          
#'     featureSets <- FeatureSetCollection(ontology)
#'     featureSets
#'     featureSets[3:5]
#'     featureSets[["Photoreceptor Activity"]]
#'     
#'     subNetworks <- list(MAPK = matrix(c("NRAS", "NRAS", "NRAS", "BRAF", "MEK",
#'                                         "ARAF", "BRAF", "CRAF", "MEK", "ERK"), ncol = 2),
#'                         P53 = matrix(c("ATM", "ATR", "ATR", "P53",
#'                                        "CHK2", "CHK1", "P53", "MDM2"), ncol = 2)
#'                         )
#'     networkSets <- FeatureSetCollection(subNetworks)                        
#'     networkSets
#'     
#' @export
#' @usage NULL
setGeneric("FeatureSetCollection", function(sets, ...)
standardGeneric("FeatureSetCollection"))

#' @usage NULL
#' @export
setMethod("FeatureSetCollection", c("list"),
          function(sets)
          {
            new("FeatureSetCollection", sets = sets)
          })

#' @usage NULL
#' @export
setMethod("length", c("FeatureSetCollection"),
    function(x)
{
    length(x@sets)
})

#' @usage NULL
#' @rdname FeatureSetCollection
#' @export
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

#' @export
setMethod("[", c("FeatureSetCollection", "numeric", "missing", "ANY"),
    function(x, i, j, ..., drop = TRUE)
{
    new("FeatureSetCollection", sets = x@sets[i])
})

#' @export
setMethod("[[", c("FeatureSetCollection", "ANY", "missing"),
    function(x, i, j, ...)
{
    x@sets[[i]]
})

setClassUnion("FeatureSetCollectionOrNULL", c("FeatureSetCollection", "NULL"))

##### SelectParams #####

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

# Parameters for Feature Selection
setGeneric("SelectParams", function(featureRanking, ...)
standardGeneric("SelectParams"))

# Default constructor.
setMethod("SelectParams", "missing", function()
{
  new("SelectParams", featureRanking = differentMeansRanking,
      characteristics = S4Vectors::DataFrame(characteristic = "Selection Name", value = "Difference in Means"),
      minPresence = 1, intermediate = character(0), subsetToSelections = TRUE,
      tuneParams = list(nFeatures = seq(10, 100, 10), performanceType = "Balanced Error"))
})

setMethod("SelectParams", c("functionOrList"),
          function(featureRanking, characteristics = DataFrame(), minPresence = 1, 
                   intermediate = character(0), subsetToSelections = TRUE, tuneParams = list(nFeatures = seq(10, 100, 10), performanceType = "Balanced Error"), ...)
          {
            if(!is.list(featureRanking) && (ncol(characteristics) == 0 || !"Selection Name" %in% characteristics[, "characteristic"]))
            {
              characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Selection Name", value = .ClassifyRenvir[["functionsTable"]][.ClassifyRenvir[["functionsTable"]][, "character"] == attr(featureRanking, "name"), "name"]))
            }
            if(is.list(featureRanking) && (ncol(characteristics) == 0 || !"Ensemble Selection" %in% characteristics[, "characteristic"]))
            {
              selectMethodNames <- unlist(lapply(featureRanking, function(rankingFunction) .ClassifyRenvir[["functionsTable"]][.ClassifyRenvir[["functionsTable"]][, "character"] == attr(rankingFunction, "name"), "name"]))
              characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Ensemble Selection", value = paste(selectMethodNames, collapse = ", ")))
            }
            others <- list(...)
            if(length(others) == 0) others <- NULL
            new("SelectParams", featureRanking = featureRanking,
                characteristics = characteristics, minPresence = minPresence,
                intermediate = intermediate, subsetToSelections = subsetToSelections,
                tuneParams = tuneParams, otherParams = others)
          })

#' Container for Storing Details of Feature Selection Function(s)
#' 
#' @rdname SelectParams-class
#' @param object An object of class \code{SelectParams} to inspect.
#' @export
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




##### TrainParams #####

setClass("TrainParams", representation(
  classifier = "function",
  characteristics = "DataFrame",
  intermediate = "character",
  tuneParams = "listOrNULL",
  otherParams = "listOrNULL",
  getFeatures = "functionOrNULL"), contains = "StageParams")

# Parameters for Classifier Training
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
              characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Classifier Name", value = .ClassifyRenvir[["functionsTable"]][.ClassifyRenvir[["functionsTable"]][, "character"] == attr(classifier, "name"), "name"]))
            }
            new("TrainParams", classifier = classifier, characteristics = characteristics,
                intermediate = intermediate, getFeatures = getFeatures, tuneParams = tuneParams,
                otherParams = list(...))
          })

#' Inspect Model Training Details
#'
#' @rdname TrainParams-class
#' @param object An object of class \code{TrainParams} to inspect.
#' @export
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
          })

##### PredictParams #####

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
            others <- list(...)
            if(length(others) == 0) others <- NULL
            new("PredictParams", predictor = predictor, characteristics = characteristics,
                intermediate = intermediate, otherParams = others)
          })

#' Inspect Prediction Function Details
#'
#' @rdname PredictParams-class
#' @param object An object of class \code{TrainParams} to inspect.
#' @export
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

setClass("ModellingParams", representation(
  balancing = "character",
  transformParams = "TransformParamsOrNULL",
  selectParams = "SelectParamsOrNULL",
  trainParams = "TrainParams",
  predictParams = "PredictParamsOrNULL",
  doImportance = "logical"
))

ModellingParams <- function(balancing = c("downsample", "upsample", "none"),
                            transformParams = NULL, selectParams = SelectParams(),
                            trainParams = TrainParams(), predictParams = PredictParams(),
                            doImportance = FALSE)
{
  balancing <- match.arg(balancing)
  new("ModellingParams", balancing = balancing, transformParams = transformParams,
      selectParams = selectParams, trainParams = trainParams, predictParams = predictParams,
      doImportance = doImportance)
}

setClassUnion("ModellingParamsOrNULL", c("ModellingParams", "NULL"))


##### ClassifyResult #####

#' Container for Storing Classification Results
#' 
#' Contains a list of models, table of actual sample classes and predicted
#' classes, the identifiers of features selected for each fold of each
#' permutation or each hold-out classification, and performance metrics such as
#' error rates. This class is not intended to be created by the user. It is
#' created by \code{\link{crossValidate}}.
#' 
#' @name ClassifyResult
#' @rdname ClassifyResult-class
#' @aliases ClassifyResult ClassifyResult-class
#' ClassifyResult,DataFrame,character,characterOrDataFrame-method
#' show,ClassifyResult-method sampleNames sampleNames,ClassifyResult-method
#' predictions predictions,ClassifyResult-method actualOutcome
#' actualOutcome,ClassifyResult-method features features,ClassifyResult-method
#' models models,ClassifyResult-method performance
#' performance,ClassifyResult-method tunedParameters
#' tunedParameters,ClassifyResult-method totalPredictions
#' totalPredictions,ClassifyResult-method
#' @docType class
#' 
#' @section Constructor:
#' \preformatted{ClassifyResult(characteristics, originalNames, originalFeatures,
#'               rankedFeatures, chosenFeatures, models, tunedParameters, predictions, actualOutcome, importance = NULL, modellingParams = NULL, finalModel = NULL)}
#' \describe{
#' \item{\code{characteristics}}{A \code{\link{DataFrame}} describing the
#' characteristics of classification done. First column must be named
#' \code{"charateristic"} and second column must be named \code{"value"}. If
#' using wrapper functions for feature selection and classifiers in this
#' package, the function names will automatically be generated and therefore it
#' is not necessary to specify them.}
#' \item{\code{originalNames}}{All sample names.}
#' \item{\code{originalFeatures}}{All feature names. Character vector
#' or \code{\link{DataFrame}} with one row for each feature if the data set has multiple kinds
#' of measurements on the same set of samples.}
#' \item{\code{chosenFeatures}}{Features selected at each fold. Character
#' vector or a data frame if data set has multiple kinds of measurements on the same set of samples.}
#' \item{\code{models}}{All of the models fitted to the training data.}
#' \item{\code{tunedParameters}}{Names of tuning parameters and the value chosen of each parameter.}
#' \item{\code{predictions}}{A data frame containing sample IDs, predicted class or risk and information about the 
#' cross-validation iteration in which the prediction was made.}
#' \item{\code{actualOutcome}}{The known class or survival data of each sample.}
#' \item{\code{importance}}{The changes in model performance for each selected variable when it is excluded.}
#' \item{\code{modellingParams}}{Stores the object used for defining the model building to enable future reuse.}
#' \item{\code{finalModel}}{A model built using all of the sample for future use. For any tuning parameters, the
#' most popular value of the parameter in cross-validation is used.}
#' }
#' 
#' @section Summary:
#' \code{result} is a \code{ClassifyResult} object.
#' \describe{
#' \item{}{
#'     \code{show(result)}: Prints a short summary of what \code{result} contains.
#' }}
#' 
#' @section Accessors:
#' \code{result} is a \code{ClassifyResult} object.
#' \describe{
#' \item{\code{sampleNames(result)}}{Returns a vector of sample names present in the data set.}}
#' \describe{
#' \item{\code{actualOutcome(result)}}{Returns the known outcome of each sample.}}
#' \describe{
#' \item{\code{models(result)}}{A \code{list} of the models fitted for each training.}}
#' \describe{
#' \item{\code{chosenFeatureNames(result)}}{A \code{list} of the features selected for each training.}}
#' \describe{
#' \item{\code{predictions(result)}}{Returns a \code{DataFrame} which has columns with test sample,
#' cross-validation and prediction information.}}
#' \describe{
#' \item{\code{performance(result)}}{Returns a \code{list} of performance measures. This is
#' empty until \code{calcCVperformance} has been used.}}
#' \describe{
#' \item{\code{tunedParameters(result)}}{Returns a \code{list} of tuned parameter values.
#' If cross-validation is used, this list will be large, as it stores chosen values
#' for every iteration.}}
#' \describe{
#' \item{\code{totalPredictions(result)}}{A single number representing the total number.
#' of predictions made during the cross-validation procedure.}}
#'        
#' @author Dario Strbenac
#' @examples
#' 
#'   #if(require(sparsediscrim))
#'   #{
#'     data(asthma)
#'     classified <- crossValidate(measurements, classes)
#'     class(classified)
#'   #}
#'   
#' @importFrom S4Vectors as.data.frame
#' @usage NULL
#' @export
setGeneric("ClassifyResult", function(characteristics, originalNames, ...)
standardGeneric("ClassifyResult"))

#' @rdname ClassifyResult-class
#' @exportClass ClassifyResult
setClass("ClassifyResult", representation(
  characteristics = "DataFrame",
  originalNames = "character",
  originalFeatures = "characterOrDataFrame",    
  rankedFeatures = "listOrNULL",
  chosenFeatures = "listOrNULL",
  actualOutcome = "factorOrSurv",
  models = "list",
  tune = "listOrNULL",
  predictions = "DataFrame",
  performance = "listOrNULL",
  importance = "DataFrameOrNULL",
  modellingParams = "ModellingParamsOrNULL",
  finalModel = "listOrNULL")
)
#' @rdname ClassifyResult-class
#' @usage NULL
#' @export
setMethod("ClassifyResult", c("DataFrame", "character"),
          function(characteristics, originalNames, originalFeatures,
                   rankedFeatures, chosenFeatures, models, tunedParameters, predictions, actualOutcome, importance = NULL, modellingParams = NULL, finalModel = NULL)
          {
            new("ClassifyResult", characteristics = characteristics, originalNames = originalNames, originalFeatures = originalFeatures,
                rankedFeatures = rankedFeatures, chosenFeatures = chosenFeatures,
                models = models, tune = tunedParameters,
                predictions = predictions, actualOutcome = actualOutcome, importance = importance, modellingParams = modellingParams, finalModel = finalModel)
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


################################################################################
#
# Set generic accessors
#
################################################################################

#' @rdname ClassifyResult-class
#' @usage NULL
#' @export
setGeneric("sampleNames", function(object, ...)
standardGeneric("sampleNames"))

#' @rdname ClassifyResult-class
#' @usage NULL
#' @export
setMethod("sampleNames", "ClassifyResult",
          function(object)
          {
            object@originalNames
          })

#' @rdname ClassifyResult-class
#' @usage NULL
#' @export
setGeneric("allFeatureNames", function(object, ...)
standardGeneric("allFeatureNames"))

#' @rdname ClassifyResult-class
#' @usage NULL
#' @export
setMethod("allFeatureNames", c("ClassifyResult"),
          function(object)
          {
            object@originalFeatures
          })

#' @rdname ClassifyResult-class
#' @usage NULL
#' @export
setGeneric("chosenFeatureNames", function(object, ...)
standardGeneric("chosenFeatureNames"))

#' @rdname ClassifyResult-class
#' @usage NULL
#' @export
setMethod("chosenFeatureNames", "ClassifyResult",
          function(object)
          {
            object@chosenFeatures
          })

#' @export
#' @usage NULL
setGeneric("models", function(object, ...)
standardGeneric("models"))

#' @rdname ClassifyResult-class
#' @usage NULL
#' @export
setMethod("models", "ClassifyResult",
          function(object)
          {
            object@models
          })

#' @export
#' @usage NULL
setGeneric("predictions", function(object, ...)
standardGeneric("predictions"))

#' @rdname ClassifyResult-class
#' @usage NULL
#' @export
setMethod("predictions", "ClassifyResult",
          function(object)
          {
            object@predictions
          })

#' @export
#' @usage NULL
setGeneric("performance", function(object, ...)
standardGeneric("performance"))

#' @rdname ClassifyResult-class
#' @usage NULL
#' @export
setMethod("performance", "ClassifyResult",
          function(object)
          {
            object@performance
          })

#' @export
#' @usage NULL
setGeneric("actualOutcome", function(object, ...)
standardGeneric("actualOutcome"))

#' @rdname ClassifyResult-class
#' @usage NULL
#' @export
setMethod("actualOutcome", "ClassifyResult",
          function(object)
          {
            object@actualOutcome
          })

#' @export
setGeneric("tunedParameters", function(object, ...)
standardGeneric("tunedParameters"))

#' @rdname ClassifyResult-class
#' @usage NULL
#' @export
setMethod("tunedParameters", "ClassifyResult",
          function(object)
          {
            object@tune
          })

#' @export
#' @usage NULL
setGeneric("totalPredictions", function(result, ...)
standardGeneric("totalPredictions"))

#' @rdname ClassifyResult-class
#' @usage NULL
#' @export
setMethod("totalPredictions", "ClassifyResult",
          function(result)
          {
              nrow(predictions(result))
          })

setClass("MixModelsListsSet", representation(set = "list"))
setGeneric("MixModelsListsSet", function(set, ...) standardGeneric("MixModelsListsSet"))
setMethod("MixModelsListsSet", c("list"), function(set) new("MixModelsListsSet", set = set))