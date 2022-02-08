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
#' Trained dlda Object
#' 
#' Enables S4 method dispatching on it.
#' 
#' 
#' @name dlda
#' @aliases dlda dlda-class
#' @docType class
#' @author Dario Strbenac
dlda <- function(x, ...) {
  UseMethod("dlda")
}

setOldClass("dlda")


#' Trained pamr Object
#' 
#' Enables S4 method dispatching on it.
#' 
#' 
#' @name pamrtrained
#' @aliases pamrtrained pamrtrained-class
#' @docType class
#' @author Dario Strbenac
setOldClass("pamrtrained")


#' Trained svm Object
#' 
#' Enables S4 method dispatching on it.
#' 
#' 
#' @name svm
#' @aliases svm svm-class
#' @docType class
#' @author Dario Strbenac
setOldClass("svm")

#' Trained multnet Object
#' 
#' Enables S4 method dispatching on it.
#' 
#' 
#' @name multnet
#' @aliases multnet multnet-class
#' @docType class
#' @author Dario Strbenac
setOldClass("multnet")



#' Trained randomForest Object
#' 
#' Enables S4 method dispatching on it.
#' 
#' 
#' @name randomForest
#' @aliases randomForest randomForest-class
#' @docType class
#' @author Dario Strbenac
setOldClass("randomForest")


#' Trained coxph Object
#' 
#' Enables S4 method dispatching on it.
#' 
#' 
#' @name randomForest
#' @aliases coxph coxph-class
#' @docType class
setOldClass("coxph")

setOldClass("Surv")

################################################################################
#
# Create union of classes
#
################################################################################





#' Union of A Function and NULL
#' 
#' Allows a slot to be either a function or empty.
#' 
#' 
#' @name functionOrNULL
#' @aliases functionOrNULL functionOrNULL-class
#' @docType class
#' @author Dario Strbenac
#' @examples
#' 
#'   PredictParams(NULL)
#'   PredictParams(DLDApredictInterface)
#' 
setClassUnion("functionOrNULL", c("function", "NULL"))

#' Union of Functions and List of Functions
#' 
#' Allows a slot to be either a function or a list of functions.
#' 
#' 
#' @name functionOrList
#' @aliases functionOrList functionOrList-class
#' @docType class
#' @author Dario Strbenac
#' @examples
#' 
#'   SelectParams(limmaRanking)
#'   SelectParams(list(limmaRanking, limmaRanking))
#' 
setClassUnion("functionOrList", c("function", "list"))




#' Union of A Numeric Value and NULL
#' 
#' Allows a slot to be either a numeric value or empty. No constructor.
#' 
#' 
#' @name numericOrNULL
#' @aliases numericOrNULL numericOrNULL-class
#' @docType class
#' @author Dario Strbenac
setClassUnion("numericOrNULL", c("numeric", "NULL"))




#' Union of a Character and a DataFrame
#' 
#' Allows a slot to be either a character or a DataFrame.
#' 
#' 
#' @name characterOrDataFrame
#' @aliases characterOrDataFrame characterOrDataFrame-class
#' @docType class
#' @author Dario Strbenac
#' @examples
#' 
#'   setClass("Selections", representation(features = "characterOrDataFrame"))
#'   selections <- new("Selections", features = c("BRAF", "NRAS"))
#'   featuresTable <- DataFrame(assay = c("RNA-seq", "Mass spectrometry"),
#'                              feature = c("PD-1", "MITF"))
#'   omicsSelections <- new("Selections", features = featuresTable)
setClassUnion("characterOrDataFrame", c("character", "DataFrame"))

setClassUnion("factorOrSurv", c("factor", "Surv"))



#' Union of a List and NULL
#' 
#' Allows a slot to be either a list or a NULL.
#' 
#' 
#' @name listOrNULL
#' @aliases listOrNULL listOrNULL-class
#' @docType class
#' @author Dario Strbenac
#' @examples
#' 
#'   setClass("EasyClassifier", representation(model = "listOrNULL"))
#'   classifier <- new("EasyClassifier", model = NULL) # Optimistic classifier.
#' 
setClassUnion("listOrNULL", c("list", "NULL"))

#' Union of A DataFrame and DataFrameList Class
#' 
#' Allows cross-validation to accept data as either a \code{DataFrame} (for a
#' single data set) or \code{DataFrameList} (for a list of tables of related
#' measurements, such as different projects measuring the same outcome and the
#' same kind of measurements). No constructor.
#' 
#' 
#' @name DataFrameOrDataFrameList
#' @aliases DataFrameOrDataFrameList DataFrameOrDataFrameList-class
#' @docType class
#' @author Dario Strbenac
setClassUnion("DataFrameOrDataFrameList", c("DataFrame", "DataFrameList"))



################################################################################
#
# Params
#
################################################################################

##### CrossValParams #####

setClass("CrossValParams", representation(
    samplesSplits = "character",
    permutations = "numericOrNULL",
    percentTest = "numericOrNULL",  
    folds = "numericOrNULL",
    leave = "numericOrNULL",
    tuneMode = "character",
    parallelParams = "BiocParallelParam"
)
)

#' Parameters for Cross-validation Specification
#' 
#' Collects and checks necessary parameters required for cross-validation by
#' \code{\link{runTests}}.
#' 
#' 
#' @name CrossValParams
#' @aliases CrossValParams CrossValParams-class
#' @docType class
#' @section Constructor: \describe{ \item{}{
#' \preformatted{CrossValParams(samplesSplits = c("Permute k-Fold", "Permute
#' Percentage Split", "Leave-k-Out", "k-Fold"), permutations = 100, percentTest
#' = 25, folds = 5, leave = 2, tuneMode = c("Resubstitution", "Nested CV"),
#' parallelParams = SerialParam())} Creates a \code{CrossValParams} object
#' based on user-specified settings with default values of different
#' classification stages geared at homoscedastic numeric data.  \describe{
#' \item{list("samplesSplits")}{Default: "Permute k-Fold". A character value
#' specifying what kind of sample splitting to do.}
#' \item{list("permutations")}{Default: 100. Number of times to permute the
#' data set before it is split into training and test sets. Only relevant if
#' \code{samplesSplits} is either \code{"Permute k-Fold"} or \code{"Permute
#' Percentage Split"}.} \item{list("percentTest")}{The percentage of the data
#' set to assign to the test set, with the remainder of the samples belonging
#' to the training set. Only relevant if \code{samplesSplits} is \code{"Permute
#' Percentage Split"}.} \item{list("folds")}{The number of approximately
#' equal-sized folds to partition the samples into. Only relevant if
#' \code{samplesSplits} is \code{"Permute k-Fold"} or \code{"k-Fold"}.}
#' \item{list("leave")}{The number of samples to generate all possible
#' combination of and use as the test set.  Only relevant if
#' \code{samplesSplits} is \code{"Leave-k-Out"}. If set to 1, it is the
#' traditional leave-one-out cross-validation, sometimes written as LOOCV.}
#' \item{list("tuneMode")}{Default: Resubstitution. The scheme to use for
#' selecting any tuning parameters.} \item{list("parallelParams")}{An instance
#' of \code{\link{BiocParallelParam}} specifying the kind of parallelisation to
#' use. Default is to use two cores less than the total number of cores the
#' computer has, if it has four or more cores, otherwise one core, as is the
#' default of \code{\link{bpparam}}. To make results fully reproducible, please
#' choose a specific back-end depending on your operating system and also set
#' \code{RNGseed} to a number.} } } }
#' @author Dario Strbenac
#' @examples
#' 
#'   CrossValParams() # Default is 100 permutations and 5 folds of each.
#'   snow <- SnowParam(workers = 4, RNGseed = 999)
#'   CrossValParams("Leave-k-Out", leave = 2, parallelParams = snow)
#'   # Fully reproducible Leave-2-out cross-validation on 4 cores,
#'   # even if feature selection or classifier use random sampling.
#' 
#' @rawNamespace exportPattern("^[[:alpha:]]+")
#' @export
CrossValParams <- function(samplesSplits = c("Permute k-Fold", "Permute Percentage Split", "Leave-k-Out", "k-Fold"),
                           permutations = 100, percentTest = 25, folds = 5, leave = 2,
                           tuneMode = c("Resubstitution", "Nested CV", "none"), parallelParams = bpparam())
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
      parallelParams = parallelParams)
}



##### StageParams #####

#' StageParams Virtual Class
#' 
#' A class for any one of \code{\link{TransformParams}},
#' \code{\link{SelectParams}}, \code{\link{TrainParams}} or
#' \code{\link{PredictParams}}. Allows a method to dispatch on any of the
#' paramter objects specifying any stage of cross-validation.
#' 
#' 
#' @name StageParams
#' @aliases StageParams StageParams-class
#' @docType class
#' @author Dario Strbenac
setClass("StageParams", representation("VIRTUAL"))


#' Union of A StageParams Object and NULL
#' 
#' Allows a slot to be either a class that has StageParams as its virtual
#' parent class or empty. No constructor.
#' 
#' 
#' @name StageParamsOrMissing
#' @aliases StageParamsOrMissing StageParamsOrMissing-class
#' @docType class
#' @author Dario Strbenac
setClassUnion("StageParamsOrMissing", c("StageParams", "missing"))


#' Union of A StageParams Object and NULL
#' 
#' Allows a slot to be either a class that has StageParams as its virtual
#' parent class or empty or NULL. No constructor.
#' 
#' 
#' @name StageParamsOrMissingOrNULL
#' @aliases StageParamsOrMissingOrNULL StageParamsOrMissingOrNULL-class
#' @docType class
#' @author Dario Strbenac
setClassUnion("StageParamsOrMissingOrNULL", c("StageParams", "missing", "NULL"))




##### TransformParams #####

setClass("TransformParams", representation(
  transform = "function",
  characteristics = "DataFrame",
  intermediate = "character",  
  otherParams = "list"), contains = "StageParams"
)

#' Union of A TransformParams Object and NULL
#' 
#' Allows a slot to be either a TransformParams class object or empty. No
#' constructor.
#' 
#' 
#' @name TransformParamsOrNULL
#' @aliases TransformParamsOrNULL TransformParamsOrNULL-class
#' @docType class
#' @author Dario Strbenac
#' @examples
#' 
#'   ModellingParams(transformParams = NULL)
#'   ModellingParams(transformParams = TransformParams(subtractFromLocation),
#'                   selectParams = SelectParams(leveneRanking))
#' 
setClassUnion("TransformParamsOrNULL", c("TransformParams", "NULL"))


#' Parameters for Data Transformation
#' 
#' Collects and checks necessary parameters required for transformation. The
#' empty constructor is for when no data transformation is desired. See
#' \code{\link{subtractFromLocation}} for an example of such a function.
#' 
#' 
#' @name TransformParams
#' @aliases TransformParams TransformParams-class TransformParams,ANY-method
#' TransformParams,function-method show,TransformParams-method
#' @docType class
#' @section Constructor: \describe{ \item{}{ \code{TransformParams(transform,
#' characteristics = DataFrame(), intermediate = character(0), ...)} Creates a
#' TransformParams object which stores the function which will do the
#' transformation and parameters that the function will use.  \describe{
#' \item{list("transform")}{A function which will do the transformation. The
#' first argument must be a \code{\link{DataFrame}} object.}
#' \item{list("characteristics")}{A \code{\link{DataFrame}} describing the
#' characteristics of data transformation to be done. First column must be
#' named \code{"charateristic"} and second column must be named \code{"value"}.
#' If using wrapper functions for data transformation in this package, the data
#' transformation name will automatically be generated and therefore it is not
#' necessary to specify it.} \item{list("intermediate")}{Character vector.
#' Names of any variables created in prior stages by \code{\link{runTest}} that
#' need to be passed to a feature selection function.} \item{list("...")}{Other
#' named parameters which will be used by the transformation function.} } } }
#' @author Dario Strbenac
#' @examples
#' 
#'   transformParams <- TransformParams(subtractFromLocation, location = "median")
#'   # Subtract all values from training set median, to obtain absolute deviations.
#' 
#' @export
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




##### FeatureSetCollection #####

setClass("FeatureSetCollection", representation(
  sets = "list")
)

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
#' @section Constructor: \describe{ \item{}{ \code{FeatureSetCollection(sets)}}
#' } \describe{ \item{list("sets")}{A named list. The names of the list
#' describe the sets and the elements of the list specify the features which
#' comprise the sets.} }
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


#' Union of a FeatureSetCollection and NULL
#' 
#' Allows a slot to be either a FeatureSetCollectionOrNULL object or empty.
#' 
#' 
#' @name FeatureSetCollectionOrNULL
#' @aliases FeatureSetCollectionOrNULL FeatureSetCollectionOrNULL-class
#' @docType class
#' @author Dario Strbenac
#' @examples
#' 
#'   TrainParams(DLDAtrainInterface, transform = NULL) # Use the input data as-is.
#' 
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


#' Union of A SelectParams Object and NULL
#' 
#' Allows a slot to be either a SelectParams class object or empty. No
#' constructor.
#' 
#' 
#' @name SelectParamsOrNULL
#' @aliases SelectParamsOrNULL SelectParamsOrNULL-class
#' @docType class
#' @author Dario Strbenac
#' @examples
#' 
#'   ModellingParams(selectParams = NULL)
#'   ModellingParams(selectParams = SelectParams(differentMeansRanking))
#' 
setClassUnion("SelectParamsOrNULL", c("SelectParams", "NULL"))


#' Parameters for Feature Selection
#' 
#' Collects and checks necessary parameters required for feature selection.
#' Either one function is specified or a list of functions to perform ensemble
#' feature selection. The empty constructor is provided for convenience.
#' 
#' 
#' @name SelectParams
#' @aliases SelectParams SelectParams-class SelectParams,missing-method
#' SelectParams,functionOrList-method
#' @docType class
#' @section Constructor: \describe{ \item{}{ \code{SelectParams()} Creates a
#' default SelectParams object. This uses either an ordinary t-test or ANOVA
#' (depending on the number of classes) and tries the top 10 to top 100
#' features in increments of 10, and picks the number of features with the best
#' resubstitution balanced error rate. Users should create an appropriate
#' \code{SelectParams} object for the characteristics of their data.  }
#' \item{}{ \preformatted{SelectParams(featureSelection, characteristics =
#' DataFrame(), minPresence = 1, intermediate = character(0),
#' subsetToSelections = TRUE, ...)} Creates a \code{SelectParams} object which
#' stores the function which will do the selection and parameters that the
#' function will use.  \describe{ \item{list("featureRanking")}{Either a
#' function which will rank the features from most promising to least promising
#' or a list of such functions. For a particular function, the first argument
#' must be an \code{\link{DataFrame}} object. The function's return value must
#' be a vector of indices.} \item{list("characteristics")}{A
#' \code{\link{DataFrame}} describing the characteristics of feature selection
#' to be done. First column must be named \code{"charateristic"} and second
#' column must be named \code{"value"}. If using wrapper functions for feature
#' selection in this package, the feature selection name will automatically be
#' generated and therefore it is not necessary to specify it.}
#' \item{list("minPresence")}{If a list of functions was provided, how many of
#' those must a feature have been selected by to be used in classification. 1
#' is equivalent to a set union and a number the same length as
#' \code{featureSelection} is equivalent to set intersection.}
#' \item{list("intermediate")}{Character vector. Names of any variables created
#' in prior stages by \code{\link{runTest}} that need to be passed to a feature
#' selection function.} \item{list("subsetToSelections")}{Whether to subset the
#' data table(s), after feature selection has been done.}
#' \item{list("...")}{Other named parameters which will be used by the
#' selection function.  If \code{featureSelection} was a list of functions,
#' this must be a list of lists, as long as \code{featureSelection}.} } } }
#' @author Dario Strbenac
#' @examples
#' 
#'   #if(require(sparsediscrim))
#'   #{
#'     SelectParams(differentMeansRanking)
#'     
#'     # Ensemble feature selection.
#'     SelectParams(list(differentMeansRanking, pairsDifferencesRanking))
#'   #}
#' 
#' @export
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




##### TrainParams #####


setClass("TrainParams", representation(
  classifier = "function",
  characteristics = "DataFrame",
  intermediate = "character",
  tuneParams = "listOrNULL",
  otherParams = "listOrNULL",
  getFeatures = "functionOrNULL"
), contains = "StageParams"
)


#' Parameters for Classifier Training
#' 
#' Collects and checks necessary parameters required for classifier training.
#' The empty constructor is provided for convenience.
#' 
#' 
#' @name TrainParams
#' @aliases TrainParams TrainParams-class TrainParams,missing-method
#' TrainParams,function-method show,TrainParams-method
#' @docType class
#' @section Constructor: \describe{ \item{}{ \code{TrainParams()} Creates a
#' default TrainParams object. The classifier function is \code{dlda} for
#' Diagonal LDA. Users should create an appropriate \code{TrainParams} object
#' for the characteristics of their data, once they are familiar with this
#' software.  } \item{}{ \preformatted{TrainParams(classifier, characteristics
#' = DataFrame(), intermediate = character(0), getFeatures = NULL, ...)}
#' Creates a TrainParams object which stores the function which will do the
#' classifier building and parameters that the function will use.  \describe{
#' \item{list("classifier")}{A function which will construct a classifier, and
#' also possibly make the predictions. The first argument must be a
#' \code{\link{DataFrame}} object. The second argument must be a vector of
#' classes. If the function also makes predictions and the value of the
#' \code{predictor} setting of \code{PredictParams} is therefore \code{NULL},
#' the third argument must be a \code{DataFrame} of test data. The function
#' must also accept a parameter named \code{verbose}. The function's return
#' value can be either a trained classifier if the function only does training
#' or a vector or data frame of class predictions if it also does prediction
#' with the test set samples.} \item{list("characteristics")}{A
#' \code{\link{DataFrame}} describing the characteristics of the classifier
#' used. First column must be named \code{"charateristic"} and second column
#' must be named \code{"value"}. If using wrapper functions for classifiers in
#' this package, a classifier name will automatically be generated and
#' therefore it is not necessary to specify it.}
#' \item{list("intermediate")}{Character vector. Names of any variables created
#' in prior stages by \code{\link{runTest}} that need to be passed to
#' \code{classifier}.} \item{list("getFeatures")}{A function may be specified
#' that extracts the selected features from the trained model. This is relevant
#' if using a classifier that does feature selection within training (e.g.
#' random forest). The function must return a list of two vectors. The first
#' vector contains the ranked features (or empty if the training algorithm
#' doesn't produce rankings) and the second vector contains the selected
#' features.} \item{list("...")}{Other named parameters which will be used by
#' the classifier.} } } }
#' @author Dario Strbenac
#' @examples
#' 
#' #if(require(sparsediscrim))
#'   trainParams <- TrainParams(DLDAtrainInterface)
#' 
#' @export
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





##### PredictParams #####

setClass("PredictParams", representation(
  predictor = "functionOrNULL",
  characteristics = "DataFrame",  
  intermediate = "character",
  otherParams = "listOrNULL"), contains = "StageParams"
)

#' Parameters for Classifier Prediction
#' 
#' Collects the function to be used for making predictions and any associated
#' parameters.
#' 
#' The function specified must return either a factor vector of class
#' predictions, or a numeric vector of scores for the second class, according
#' to the levels of the class vector of the input data set, or a data frame
#' which has two columns named class and score.
#' 
#' 
#' @name PredictParams
#' @aliases PredictParams PredictParams-class PredictParams,missing-method
#' PredictParams,functionOrNULL-method show,PredictParams-method
#' @docType class
#' @section Constructor: \describe{ \item{}{ \code{PredictParams()} Creates a
#' default PredictParams object. This assumes that the object returned by the
#' classifier has a list element named \code{"class"}.  } \item{}{
#' \code{PredictParams(predictor, characteristics = DataFrame(), intermediate =
#' character(0), ...)} Creates a PredictParams object which stores the function
#' which will do the class prediction, if required, and parameters that the
#' function will use. If the training function also makes predictions, this
#' must be set to \code{NULL}.  } \describe{ \item{list("predictor")}{Either
#' \code{NULL} or a function to make predictions with. If it is a function,
#' then the first argument must accept the classifier made in the training
#' step.  The second argument must accept a \code{\link{DataFrame}} of new
#' data.} \item{list("characteristics")}{A \code{\link{DataFrame}} describing
#' the characteristics of the predictor function used. First column must be
#' named \code{"charateristic"} and second column must be named
#' \code{"value"}.} \item{list("intermediate")}{Character vector. Names of any
#' variables created in prior stages in \code{\link{runTest}} that need to be
#' passed to the prediction function.} \item{list("...")}{Other arguments that
#' \code{predictor} may use.} } }
#' @author Dario Strbenac
#' @examples
#' 
#' predictParams <- PredictParams(predictor = DLDApredictInterface)
#' # For prediction by trained object created by DLDA training function.
#' PredictParams(predictor = NULL)
#' # For when the training function also does prediction and directly returns the
#' # predictions.
#' 
#' @export
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






#' Union of A PredictParams Object and NULL
#' 
#' Allows a slot to be either a PredictParams class object or empty. No
#' constructor. In other words, the training function specified also makes
#' predictions with the test set.
#' 
#' 
#' @name PredictParamsOrNULL
#' @aliases PredictParamsOrNULL PredictParamsOrNULL-class
#' @docType class
#' @author Dario Strbenac
#' @examples
#' 
#'   ModellingParams(trainParams = TrainParams(kNNinterface, k = 5), predictParams = NULL)
#' 
setClassUnion("PredictParamsOrNULL", c("PredictParams", "NULL"))

setClass("ModellingParams", representation(
  balancing = "character",
  transformParams = "TransformParamsOrNULL",
  selectParams = "SelectParamsOrNULL",
  trainParams = "TrainParams",
  predictParams = "PredictParamsOrNULL"
))





##### ModellingParams #####

#' Parameters for Data Modelling Specification
#' 
#' Collects and checks necessary parameters required for data modelling. Apart
#' from data transfomation that needs to be done within cross-validation (e.g.
#' \code{\link{subtractFromLocation}}), feature selection, model training and
#' prediction, this container also stores a setting for class imbalance
#' rebalancing.
#' 
#' 
#' @name ModellingParams
#' @aliases ModellingParams ModellingParams-class
#' @docType class
#' @section Constructor: \describe{ \item{}{ \code{ModellingParams()} Creates a
#' default \code{ModellingParams} object suitable for classification of
#' homoscedastic numeric data.  This uses either an ordinary t-test or ANOVA
#' (depending on the number of classes) and tries the top 10 to top 100
#' features in increments of 10 during feature selection and diagonal LDA as
#' the classifier. Also, the default class rebalancing scheme is downsampling
#' to the smallest class. Users should create an appropriate
#' \code{ModellingParams} object for the characteristics of their data.  }
#' \item{}{ \preformatted{ModellingParams(balancing = c("downsample",
#' "upsample", "none"), transformParams = NULL, selectParams = SelectParams(),
#' trainParams = TrainParams(), predictParams = PredictParams())} Creates a
#' \code{ModellingParams} object based on user-specified settings with default
#' values of different classification stages geared at homoscedastic numeric
#' data.  \describe{ \item{list("balancing")}{Default: "downsample". A
#' character value specifying what kind of class balancing to do, if any.}
#' \item{list("transformParams")}{Parameters used for feature transformation
#' specified by a \code{\link{TransformParams}} instance.  Optional, can be
#' \code{NULL}}.  \item{list("selectParams")}{Parameters used during feature
#' selection specified by a \code{\link{SelectParams}} instance.  By default,
#' parameters for selection based on differences in means of numeric data.
#' Optional, can be \code{NULL}.} \item{list("trainParams")}{Parameters for
#' model training specified by a \code{\link{TrainParams}} instance.  By
#' default, uses diagonal LDA.} \item{list("predictParams")}{Parameters for
#' model training specified by a \code{\link{PredictParams}} instance.  By
#' default, uses diagonal LDA.} } } }
#' @author Dario Strbenac
#' @examples
#' 
#'   #if(require(sparsediscrim))
#'   #{
#'      ModellingParams() # Default is differences in means selection and DLDA.
#'      ModellingParams(selectParams = NULL, # No feature selection before training.
#'                      trainParams = TrainParams(randomForestTrainInterface),
#'                      predictParams = PredictParams(randomForestPredictInterface))
#'   #}
#' @export
ModellingParams <- function(balancing = c("downsample", "upsample", "none"),
                            transformParams = NULL, selectParams = SelectParams(),
                            trainParams = TrainParams(), predictParams = PredictParams())
{
  balancing <- match.arg(balancing)
  new("ModellingParams", balancing = balancing, transformParams = transformParams,
      selectParams = selectParams, trainParams = trainParams, predictParams = predictParams)
}


#' Union of A ModellingParams Object and NULL
#' 
#' Allows a slot to be either a ModellingParams class object or empty. No
#' constructor.
#' 
#' 
#' @name ModellingParamsOrNULL
#' @aliases ModellingParamsOrNULL ModellingParamsOrNULL-class
#' @docType class
setClassUnion("ModellingParamsOrNULL", c("ModellingParams", "NULL"))


##### ClassifyResult #####

#' Container for Storing Classification Results
#' 
#' Contains a list of models, table of actual sample classes and predicted
#' classes, the identifiers of features selected for each fold of each
#' permutation or each hold-out classification, and performance metrics such as
#' error rates.  This class is not intended to be created by the user. It is
#' created by \code{\link{runTest}} or \code{\link{runTests}}.
#' 
#' 
#' @name ClassifyResult
#' @aliases ClassifyResult ClassifyResult-class
#' ClassifyResult,DataFrame,character,characterOrDataFrame-method
#' show,ClassifyResult-method sampleNames sampleNames,ClassifyResult-method
#' featureNames featureNames,ClassifyResult-method predictions
#' predictions,ClassifyResult-method actualClasses
#' actualClasses,ClassifyResult-method features features,ClassifyResult-method
#' models models,ClassifyResult-method performance
#' performance,ClassifyResult-method tunedParameters
#' tunedParameters,ClassifyResult-method totalPredictions
#' totalPredictions,ClassifyResult-method
#' @docType class
#' @section Constructor: \describe{ \item{}{
#' \preformatted{ClassifyResult(characteristics, originalNames,
#' originalFeatures, rankedFeatures, chosenFeatures, predictions,
#' actualClasses, models, validation, tune = NULL)} } } \describe{
#' \item{list("characteristics")}{A \code{\link{DataFrame}} describing the
#' characteristics of classification done. First column must be named
#' \code{"charateristic"} and second column must be named \code{"value"}. If
#' using wrapper functions for feature selection and classifiers in this
#' package, the function names will automatically be generated and therefore it
#' is not necessary to specify them.} \item{list("originalNames")}{All sample
#' names.} \item{list("originalFeatures")}{All feature names. Character vector
#' or \code{\link{DataFrame}} with one row for each feature if the data set is
#' a \code{\link{MultiAssayExperiment}}.} \item{list("rankedFeatures")}{All
#' features, from most to least important. Character vector or
#' \code{\link{DataFrame}} if data set is a \code{MultiAssayExperiment}.}
#' \item{list("chosenFeatures")}{Features selected at each fold. Character
#' vector or DataFrame if data set is a \code{MultiAssayExperiment}.}
#' \item{list("predictions")}{A \code{\link{list}} of \code{\link{data.frame}}
#' containing information about samples, their actual class and predicted class
#' and/or class score and information about the cross-validation fold in which
#' the prediction was made.} \item{list("actualClasses")}{Factor of class of
#' each sample.} \item{list("models")}{All of the models fitted to the training
#' data.} \item{list("validation")}{List with first element being the name of
#' the validation scheme, and other elements providing details about the
#' scheme.} \item{list("tune")}{A description of the tuning parameters, and the
#' value chosen of each parameter.} }
#' @author Dario Strbenac
#' @examples
#' 
#'   #if(require(sparsediscrim))
#'   #{
#'     data(asthma)
#'     
#'     LOOCVparams <- CrossValParams("Leave-k-Out", leave = 1)
#'     modellingParams <- ModellingParams()
#'     classified <-
#'     runTests(measurements, classes, LOOCVparams, modellingParams,
#'              DataFrame(characteristic = c("dataset", "classification"),
#'                       value = c("Asthma", "Different Means"))
#'              )
#'     class(classified)
#'   #}
#'   
#' @importFrom S4Vectors as.data.frame  
#' @export
setGeneric("ClassifyResult", function(characteristics, originalNames, originalFeatures, ...)
standardGeneric("ClassifyResult"))
setClass("ClassifyResult", representation(
  characteristics = "DataFrame",
  originalNames = "character",
  originalFeatures = "characterOrDataFrame",
  rankedFeatures = "listOrNULL",
  chosenFeatures = "listOrNULL",
  actualClasses = "factorOrSurv",
  models = "list",
  tune = "listOrNULL",
  predictions = "data.frame",
  performance = "listOrNULL",
  modellingParams = "ModellingParamsOrNULL",
  finalModel = "listOrNULL")
)
setMethod("ClassifyResult", c("DataFrame", "character", "characterOrDataFrame"),
          function(characteristics, originalNames, originalFeatures,
                   rankedFeatures, chosenFeatures, models, tunedParameters, predictions, actualClasses, modellingParams = NULL, finalModel = NULL)
          {
            new("ClassifyResult", characteristics = characteristics,
                originalNames = originalNames, originalFeatures = originalFeatures,
                rankedFeatures = rankedFeatures, chosenFeatures = chosenFeatures,
                models = models, tune = tunedParameters,
                predictions = predictions, actualClasses = actualClasses, modellingParams = modellingParams, finalModel = finalModel)
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

#' @export
setGeneric("sampleNames", function(object, ...)
standardGeneric("sampleNames"))
setMethod("sampleNames", c("ClassifyResult"),
          function(object)
          {
            object@originalNames
          })

#' @export
setGeneric("featureNames", function(object, ...)
standardGeneric("featureNames"))
setMethod("featureNames", c("ClassifyResult"),
          function(object)
          {
            object@originalFeatures
          })

#' @export
setGeneric("features", function(object, ...)
standardGeneric("features"))
setMethod("features", c("ClassifyResult"),
          function(object)
          {
            object@chosenFeatures
          })

#' @export
setGeneric("models", function(object, ...)
standardGeneric("models"))
setMethod("models", c("ClassifyResult"),
          function(object)
          {
            object@models
          })

#' @export
setGeneric("predictions", function(object, ...)
standardGeneric("predictions"))
setMethod("predictions", c("ClassifyResult"),
          function(object)
          {
            object@predictions
          })

#' @export
setGeneric("performance", function(object, ...)
standardGeneric("performance"))
setMethod("performance", c("ClassifyResult"),
          function(object)
          {
            object@performance
          })

#' @export
setGeneric("actualClasses", function(object, ...)
standardGeneric("actualClasses"))
setMethod("actualClasses", c("ClassifyResult"),
          function(object)
          {
            object@actualClasses
          })

#' @export
setGeneric("tunedParameters", function(object, ...)
standardGeneric("tunedParameters"))
setMethod("tunedParameters", "ClassifyResult",
          function(object)
          {
            object@tune
          })

#' @export
setGeneric("totalPredictions", function(result, ...)
standardGeneric("totalPredictions"))
setMethod("totalPredictions", c("ClassifyResult"),
          function(result)
          {
              nrow(predictions(result))
          })


#' 
#' setClass("MixModelsListsSet", representation(
#'   set = "list")
#' )
#' 
#' #' Container for a List of Lists Containing Mixture Models
#' #' 
#' #' Stores a list of lists of trained mixture models, to prevent them being
#' #' unintentionally being unlisted during cross-validation. Not intended for
#' #' end-user.
#' #' 
#' #' 
#' #' @name MixModelsListsSet
#' #' @aliases MixModelsListsSet MixModelsListsSet-class
#' #' MixModelsListsSet,list-method
#' #' @docType class
#' #' @section Constructor: \describe{ \item{}{
#' #' \preformatted{MixModelsListsSet(set)} Creates a MixModelsListsSet object
#' #' which stores the mixture models.  \describe{ \item{list("set")}{A list as
#' #' long as the number of classes in the data set. Each element is a list, which
#' #' each element of is a mixture model trained on one feature.} } } }
#' #' @author Dario Strbenac
#' #' @examples
#' #' 
#' #'   if(require(Rmixmod))
#' #'   {
#' #'     mixModels <- list(Good = list(mixmodCluster(rnorm(20), nbCluster = 1:2)),
#' #'                       Poor = list(mixmodCluster(rnorm(20), nbCluster = 1:2)))
#' #'     MixModelsListsSet(mixModels)
#' #'   }
#' #' 
#' #' @export
#' setGeneric("MixModelsListsSet", function(set, ...)
#' standardGeneric("MixModelsListsSet"))
#' setMethod("MixModelsListsSet", c("list"),
#'           function(set)
#'           {
#'             new("MixModelsListsSet", set = set)
#'           })
#' 
