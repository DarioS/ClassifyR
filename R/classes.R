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
#' @name dlda-class
#' @aliases dlda dlda-class
#' @docType class
#' @author Dario Strbenac
#' @usage NULL
dlda <- function(x, ...) {
  UseMethod("dlda")
}

setOldClass("dlda")


#' Trained pamr Object
#' 
#' Enables S4 method dispatching on it.
#' 
#' 
#' @name pamrtrained-class
#' @aliases pamrtrained pamrtrained-class
#' @docType class
#' @author Dario Strbenac
setOldClass("pamrtrained")


#' Trained svm Object
#' 
#' Enables S4 method dispatching on it.
#' 
#' 
#' @name svm-class
#' @aliases svm svm-class
#' @docType class
#' @author Dario Strbenac
setOldClass("svm")

#' Trained multnet Object
#' 
#' Enables S4 method dispatching on it.
#' 
#' 
#' @name multnet-class
#' @aliases multnet multnet-class
#' @docType class
#' @author Dario Strbenac
setOldClass("multnet")

#' @name coxnet-class
#' @aliases coxnet coxnet-class
#' @docType class
#' @author Dario Strbenac
setOldClass("coxnet")



#' Trained randomForest Object
#' 
#' Enables S4 method dispatching on it.
#' 
#' 
#' @name randomForest-class
#' @aliases randomForest randomForest-class
#' @docType class
#' @author Dario Strbenac
setOldClass("randomForest")


#' Trained coxph Object
#' 
#' Enables S4 method dispatching on it.
#' 
#' 
#' @name coxph-class
#' @aliases coxph coxph-class
#' @docType class
setOldClass("coxph")

#' Survival Data Container
#' 
#' Enables S4 method dispatching on it.
#' 
#' 
#' @name Surv-class
#' @aliases Surv Surv-class
#' @docType class
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
#' @name functionOrNULL-class
#' @aliases functionOrNULL functionOrNULL-class
#' @docType class
#' @author Dario Strbenac
#' @examples
#' 
#'   PredictParams(NULL) # Training function does both tasks.
#'   PredictParams(DLDApredictInterface)
#' 
setClassUnion("functionOrNULL", c("function", "NULL"))

#' Union of Functions and List of Functions
#' 
#' Allows a slot to be either a function or a list of functions.
#' 
#' 
#' @name functionOrList-class
#' @aliases functionOrList functionOrList-class
#' @docType class
#' @author Dario Strbenac
#' @examples
#' 
#'   SelectParams(limmaRanking)
#'   SelectParams(list(limmaRanking, differentMeansRanking)) # Ensemble selection.
#' 
setClassUnion("functionOrList", c("function", "list"))


#' Union of A Numeric Value and NULL
#' 
#' Allows a slot to be either a numeric value or empty. No constructor.
#' 
#' 
#' @name numericOrNULL-class
#' @aliases numericOrNULL numericOrNULL-class
#' @docType class
#' @author Dario Strbenac
setClassUnion("numericOrNULL", c("numeric", "NULL"))


#' Union of a Character and a DataFrame
#' 
#' Allows a slot to be either a character or a DataFrame.
#' 
#' 
#' @name characterOrDataFrame-class
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

#' Union of a Surv class and a factor
#' 
#' Enables different functionality to be executed depending on whether
#' input data dependent variable is survival or categorical classes.
#' 
#' 
#' @name characterOrDataFrame-class
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
setClassUnion("factorOrSurv", c("factor", "Surv"))



#' Union of a List and NULL
#' 
#' Allows a slot to be either a list or a NULL.
#' 
#' 
#' @name listOrNULL-class
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
#' @name DataFrameOrDataFrameList-class
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

#' Parameters for Cross-validation Specification
#' 
#' Collects and checks necessary parameters required for cross-validation by
#' \code{\link{runTests}}.
#' 
#' 
#' @name CrossValParams
#' @rdname CrossValParams-class
#' @aliases CrossValParams CrossValParams-class
#' @docType class
#' 
#' @param samplesSplits Default: "Permute k-Fold". A character value
#' specifying what kind of sample splitting to do.
#' @param permutations Default: 100. Number of times to permute the
#' data set before it is split into training and test sets. Only relevant if
#' \code{samplesSplits} is either \code{"Permute k-Fold"} or \code{"Permute
#' Percentage Split"}.
#' @param percentTest The percentage of the data
#' set to assign to the test set, with the remainder of the samples belonging
#' to the training set. Only relevant if \code{samplesSplits} is \code{"Permute
#' Percentage Split"}.
#' @param folds The number of approximately equal-sized folds to partition
#' the samples into. Only relevant if \code{samplesSplits} is \code{"Permute k-Fold"}
#' or \code{"k-Fold"}.
#' @param leave The number of samples to generate all possible
#' combination of and use as the test set.  Only relevant if \code{samplesSplits} is
#' \code{"Leave-k-Out"}. If set to 1, it is the traditional leave-one-out cross-validation,
#' sometimes written as LOOCV.
#' @param tuneMode Default: Resubstitution. The scheme to use for selecting any tuning parameters.
#' @param parallelParams An instance of \code{\link{BiocParallelParam}} specifying
#' the kind of parallelisation to use. Default is to use two cores less than the total number of
#' cores the computer has, if it has four or more cores, otherwise one core, as is the
#' default of \code{\link{bpparam}}. To make results fully reproducible, please
#' choose a specific back-end depending on your operating system and also set
#' \code{RNGseed} to a number.
#' 
#' @author Dario Strbenac
#' @examples
#' 
#'   CrossValParams() # Default is 100 permutations and 5 folds of each.
#'   snow <- SnowParam(workers = 4, RNGseed = 999)
#'   CrossValParams("Leave-k-Out", leave = 2, parallelParams = snow)
#'   # Fully reproducible Leave-2-out cross-validation on 4 cores,
#'   # even if feature selection or classifier use random sampling.
#' 
#' @exportClass CrossValParams
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

# CrossValParams constructor is an ordinary function for performance reasons.
#' @export
#' @rdname CrossValParams-class
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
#' parameter objects specifying any stage of cross-validation.
#' 
#' 
#' @name StageParams-class
#' @aliases StageParams StageParams-class
#' @author Dario Strbenac
setClass("StageParams", representation("VIRTUAL"))


#' Union of A StageParams Object and NULL
#' @name StageParamsOrMissing-class
#' @rdname StageParamsOrElse
#' @aliases StageParamsOrMissing StageParamsOrMissing-class
#' @author Dario Strbenac
setClassUnion("StageParamsOrMissing", c("StageParams", "missing"))


#' Union of A StageParams Object and NULL
#' 
#' StageParamsOrMissing: Allows a slot to be either a class that has StageParams as its virtual
#' parent class or empty. No constructor.
#' StageParamsOrMissingOrNULL: Allows a slot to be either a class that has StageParams as its virtual
#' parent class or empty or NULL. No constructor.
#' 
#' @name StageParamsOrMissingOrNULL-class
#' @rdname StageParamsOrElse
#' @aliases StageParamsOrMissingOrNULL StageParamsOrMissingOrNULL-class
setClassUnion("StageParamsOrMissingOrNULL", c("StageParams", "missing", "NULL"))


##### TransformParams #####
#' @exportClass TransformParams
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
#' @name TransformParamsOrNULL-class
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
#' Collects and checks necessary parameters required for transformation within CV. The
#' empty constructor is for when no data transformation is desired. See
#' \code{\link{subtractFromLocation}} for an example of such a function.
#' 
#' 
#' @name TransformParams
#' @rdname TransformParams-class
#' @aliases TransformParams TransformParams-class TransformParams,ANY-method
#' TransformParams,function-method show,TransformParams-method
#' @docType class
#' @usage NULL
#' @section Constructor:
#' \describe{
#' \item{}{
#' \code{TransformParams(transform, characteristics = DataFrame(), intermediate = character(0), ...)} 
#' Creates a \code{TransformParams} object which stores the function which will do the
#' transformation and parameters that the function will use.
#' \describe{
#' \item{\code{transform}}{A function which will do the transformation. The
#' first argument must be a \code{\link{DataFrame}} object.}
#' \item{\code{characteristics}}{A \code{\link{DataFrame}} describing the
#' characteristics of data transformation to be done. First column must be
#' named \code{"charateristic"} and second column must be named \code{"value"}.
#' If using wrapper functions for data transformation in this package, the data
#' transformation name will automatically be generated and therefore it is not
#' necessary to specify it.}
#' \item{\code{intermediate}}{Character vector. Names of any variables created in
#' prior stages by \code{\link{runTest}} that need to be passed to a feature selection
#' function.}
#' \item{\code{...}}{Other named parameters which will be used by the transformation function.}
#' } } }
#' 
#' @section Summary:
#' \code{transformParams} is a \code{TransformParams} object.
#' \describe{
#' \item{}{
#'     \code{show(transformParams)}: Prints a short summary of what \code{transformParams} contains.
#'  }}
#' 
#' @author Dario Strbenac
#' @examples
#' 
#'   transformParams <- TransformParams(subtractFromLocation, location = "median")
#'   # Subtract all values from training set median, to obtain absolute deviations.
#' 
#' @export
#' @usage NULL
setGeneric("TransformParams", function(transform, ...)
standardGeneric("TransformParams"))

#' @rdname TransformParams-class
#' @usage NULL
#' @export
setMethod("TransformParams", "function",
          function(transform, characteristics = S4Vectors::DataFrame(), intermediate = character(0), ...)
          {
            if(ncol(characteristics) == 0 || !"Transform Name" %in% characteristics[, "characteristic"])
            {
              characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Transform Name", value = .ClassifyRenvir[["functionsTable"]][.ClassifyRenvir[["functionsTable"]][, "character"] == transform@generic, "name"]))
            }
            new("TransformParams", transform = transform, characteristics = characteristics,
                intermediate = intermediate, otherParams = list(...))
          })

#' @usage NULL
#' @rdname TransformParams-class
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

#' @usage NULL
#' @export
setMethod("[", c("FeatureSetCollection", "numeric", "missing", "ANY"),
    function(x, i, j, ..., drop = TRUE)
{
    new("FeatureSetCollection", sets = x@sets[i])
})

#' @usage NULL
#' @export
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
#' @name FeatureSetCollectionOrNULL-class
#' @aliases FeatureSetCollectionOrNULL FeatureSetCollectionOrNULL-class
#' @docType class
#' @author Dario Strbenac
#' @examples
#' 
#'   TrainParams(DLDAtrainInterface, transform = NULL) # Use the input data as-is.
#' 
setClassUnion("FeatureSetCollectionOrNULL", c("FeatureSetCollection", "NULL"))





##### SelectParams #####

#' @exportClass SelectParams
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
#' @name SelectParamsOrNULL-class
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
#' @rdname SelectParams-class
#' @aliases SelectParams SelectParams-class SelectParams,missing-method
#' SelectParams,functionOrList-method
#' @docType class
#' @section Constructor:
#' \describe{ \item{}{
#' \code{SelectParams()} Creates a default \code{SelectParams} object. This uses either
#' an ordinary t-test or ANOVA (depending on the number of classes) and tries the
#' top 10 to top 100 features in increments of 10, and picks the number of features
#' with the best resubstitution balanced error rate. Users should create an appropriate
#' \code{SelectParams} object for the characteristics of their data.}
#' \item{}{\preformatted{SelectParams(featureSelection, characteristics = DataFrame(), minPresence = 1, intermediate = character(0),
#' subsetToSelections = TRUE, tuneParams = list(nFeatures = seq(10, 100, 10), performanceType = "Balanced Error"), ...)} Creates a \code{SelectParams}
#' object which stores the function(s) which will do the selection and parameters that the
#' function will use.
#' \describe{\item{\code{featureRanking}}{Either a function which will rank the features
#' from most promising to least promising or a list of such functions. For a particular function,
#' the first argument must be an \code{\link{DataFrame}} object. The function's return value must
#' be a vector of indices.}
#' \item{\code{characteristics}}{A \code{\link{DataFrame}} describing the characteristics
#' of feature selection to be done. First column must be named \code{"charateristic"} and
#' second column must be named \code{"value"}. If using wrapper functions for feature
#' selection in this package, the feature selection name will automatically be
#' generated and therefore it is not necessary to specify it.}
#' \item{\code{minPresence}}{If a list of functions was provided, how many of
#' those must a feature have been selected by to be used in classification. 1
#' is equivalent to a set union and a number the same length as
#' \code{featureSelection} is equivalent to set intersection.}
#' \item{\code{intermediate}}{Character vector. Names of any variables created
#' in prior stages by \code{\link{runTest}} that need to be passed to a feature
#' selection function.}
#' \item{\code{subsetToSelections}}{Whether to subset the data table(s), after feature selection has been done.}
#' \item{\code{tuneParams}}{A list specifying tuning parameters required during feature selection. The names of
#' the list are the names of the parameters and the vectors are the values of the parameters to try. All possible
#' combinations are generated. Two elements named \code{nFeatures} and \code{performanceType} are mandatory, to
#' define the performance metric which will be used to select features and how many top-ranked features to try.}
#' \item{\code{...}}{Other named parameters which will be used by the
#' selection function. If \code{featureSelection} was a list of functions,
#' this must be a list of lists, as long as \code{featureSelection}.} } } }
#' @section Summary:
#' \code{selectParams} is a \code{SelectParams} object.
#' \describe{
#' \item{}{
#'   \code{show(SelectParams)}: Prints a short summary of what \code{selectParams} contains.
#' }}
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
#' @usage NULL
setGeneric("SelectParams", function(featureRanking, ...)
standardGeneric("SelectParams"))

#' @rdname SelectParams-class
#' @usage NULL
#' @export
setMethod("SelectParams", "missing", function()
{
  new("SelectParams", featureRanking = differentMeansRanking,
      characteristics = S4Vectors::DataFrame(characteristic = "Selection Name", value = "Difference in Means"),
      minPresence = 1, intermediate = character(0), subsetToSelections = TRUE,
      tuneParams = list(nFeatures = seq(10, 100, 10), performanceType = "Balanced Error"))
})
#' @rdname SelectParams-class
#' @usage NULL
#' @export
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

#' @usage NULL
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

#' @exportClass TrainParams
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
#' @rdname TrainParams-class
#' @aliases TrainParams TrainParams-class TrainParams,missing-method
#' TrainParams,function-method show,TrainParams-method
#' @docType class
#' @section Constructor:
#' \describe{
#' \item{}{\code{TrainParams()} Creates a default \code{TrainParams} object.
#' The classifier function is \code{dlda} for Diagonal LDA. Users should create
#' an appropriate \code{TrainParams} object for the characteristics of their data,
#' once they are familiar with this software.}
#' \item{}{\preformatted{TrainParams(classifier, characteristics = DataFrame(),
#' intermediate = character(0), getFeatures = NULL, ...)}
#' Creates a \code{TrainParams} object which stores the function which will do the
#' classifier building and parameters that the function will use.
#' \describe{
#' \item{\code{classifier}}{A function which will construct a classifier, and
#' also possibly make the predictions. The first argument must be a
#' \code{\link{DataFrame}} object. The second argument must be a vector of
#' classes. If the function also makes predictions and the value of the
#' \code{predictor} setting of \code{PredictParams} is therefore \code{NULL},
#' the third argument must be a \code{DataFrame} of test data. The function
#' must also accept a parameter named \code{verbose}. The function's return
#' value can be either a trained classifier if the function only does training
#' or a vector or data frame of class predictions if it also does prediction
#' with the test set samples.}
#' \item{\code{characteristics}}{A \code{\link{DataFrame}} describing the
#' characteristics of the classifier used. First column must be named \code{"charateristic"}
#' and second column must be named \code{"value"}. If using wrapper functions for classifiers
#' in this package, a classifier name will automatically be generated and
#' therefore it is not necessary to specify it.}
#' \item{\code{intermediate}}{Character vector. Names of any variables created
#' in prior stages by \code{\link{runTest}} that need to be passed to
#' \code{classifier}.}
#' \item{\code{getFeatures}}{A function may be specified that extracts the selected
#' features from the trained model. This is relevant if using a classifier that does
#' feature selection within training (e.g. random forest). The function must return a
#' list of two vectors. The first vector contains the ranked features (or empty if the
#' training algorithm doesn't produce rankings) and the second vector contains the selected
#' features.}
#' \item{\code{...}}{Other named parameters which will be used by the classifier.} } } }
#' @section Summary:
#' \code{trainParams} is a \code{TrainParams} object.
#' \describe{
#' \item{}{
#'   \code{show(trainParams)}: Prints a short summary of what \code{trainParams} contains.
#' }}
#' @author Dario Strbenac
#' @examples
#' 
#' #if(require(sparsediscrim))
#'   trainParams <- TrainParams(DLDAtrainInterface)
#' 
#' @usage NULL
#' @export
setGeneric("TrainParams", function(classifier, ...) standardGeneric("TrainParams"))

#' @usage NULL
#' @rdname TrainParams-class
#' @export
setMethod("TrainParams", "missing", function()
{
  new("TrainParams", classifier = DLDAtrainInterface,
      characteristics = S4Vectors::DataFrame(characteristic = "Classifier Name", value = "Diagonal LDA"),
      intermediate = character(0), getFeatures = NULL)
})

#' @usage NULL
#' @rdname TrainParams-class
#' @export
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

#' @usage NULL
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

            if(!is.null(object@getFeatures))
              cat("Selected Features Extracted By: ", object@getFeatures@generic, ".\n", sep = '')
          })





##### PredictParams #####

#' @exportClass PredictParams
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
#' @rdname PredictParams-class
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
#' must be set to \code{NULL}.  } \describe{ \item{\code{predictor}}{Either
#' \code{NULL} or a function to make predictions with. If it is a function,
#' then the first argument must accept the classifier made in the training
#' step.  The second argument must accept a \code{\link{DataFrame}} of new
#' data.} \item{\code{characteristics}}{A \code{\link{DataFrame}} describing
#' the characteristics of the predictor function used. First column must be
#' named \code{"charateristic"} and second column must be named
#' \code{"value"}.} \item{\code{intermediate}}{Character vector. Names of any
#' variables created in prior stages in \code{\link{runTest}} that need to be
#' passed to the prediction function.} \item{\code{...}}{Other arguments that
#' \code{predictor} may use.} } }
#' @section Summary:
#' \code{predictParams} is a \code{PredictParams} object.
#' \describe{
#' \item{}{
#'   \code{show(predictParams)}: Prints a short summary of what \code{predictParams} contains.
#' }}
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
#' @usage NULL
setGeneric("PredictParams", function(predictor, ...)
standardGeneric("PredictParams"))

#' @usage NULL
#' @rdname PredictParams-class
#' @export
setMethod("PredictParams", "missing", function()
{
  new("PredictParams", predictor = DLDApredictInterface,
      characteristics = S4Vectors::DataFrame(characteristic = "Predictor Name", value = "Diagonal LDA"),
      intermediate = character(0), otherParams = NULL)
})

#' @usage NULL
#' @rdname PredictParams-class
#' @export
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

#' @usage NULL
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

#' @exportClass ModellingParams
setClass("ModellingParams", representation(
  balancing = "character",
  transformParams = "TransformParamsOrNULL",
  selectParams = "SelectParamsOrNULL",
  trainParams = "TrainParams",
  predictParams = "PredictParamsOrNULL",
  doImportance = "logical"
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
#' @name ModellingParams
#' @rdname ModellingParams-class
#' @aliases ModellingParams ModellingParams-class
#' @docType class
#' @param balancing Default: "downsample". A character value specifying what kind
#' of class balancing to do, if any.
#' @param transformParams Parameters used for feature transformation inside of C.V.
#' specified by a \code{\link{TransformParams}} instance. Optional, can be \code{NULL}.
#' @param selectParams Parameters used during feature selection specified
#'   by a \code{\link{SelectParams}} instance.  By default, parameters for selection
#'   based on differences in means of numeric data. Optional, can be \code{NULL}.
#' @param trainParams Parameters for model training specified by a \code{\link{TrainParams}} instance.
#'   By default, uses diagonal LDA.
#' @param predictParams Parameters for model training specified by a \code{\link{PredictParams}} instance.
#' By default, uses diagonal LDA.
#' @param doImportance Default: \code{TRUE}. Whether or not to carry out removal of each feature, one at a time, which
#' was chosen and then retrain and model and predict the test set, to measure the change in performance metric. Can
#' also be set to FALSE if not of interest to reduce the modelling run time.
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
                            trainParams = TrainParams(), predictParams = PredictParams(),
                            doImportance = TRUE)
{
  balancing <- match.arg(balancing)
  new("ModellingParams", balancing = balancing, transformParams = transformParams,
      selectParams = selectParams, trainParams = trainParams, predictParams = predictParams,
      doImportance = doImportance)
}


#' Union of A ModellingParams Object and NULL
#' 
#' Allows a slot to be either a ModellingParams class object or empty. No
#' constructor.
#' 
#' 
#' @name ModellingParamsOrNULL-class
#' @aliases ModellingParamsOrNULL ModellingParamsOrNULL-class
#' @docType class
setClassUnion("ModellingParamsOrNULL", c("ModellingParams", "NULL"))


##### ClassifyResult #####

#' Container for Storing Classification Results
#' 
#' Contains a list of models, table of actual sample classes and predicted
#' classes, the identifiers of features selected for each fold of each
#' permutation or each hold-out classification, and performance metrics such as
#' error rates. This class is not intended to be created by the user. It is
#' created by \code{\link{runTest}} or \code{\link{runTests}}.
#' 
#' @name ClassifyResult
#' @rdname ClassifyResult-class
#' @aliases ClassifyResult ClassifyResult-class
#' ClassifyResult,DataFrame,character,characterOrDataFrame-method
#' show,ClassifyResult-method sampleNames sampleNames,ClassifyResult-method
#' allFeatureNames allFeatureNames,ClassifyResult-method
#' predictions predictions,ClassifyResult-method actualOutcomes
#' actualOutcomes,ClassifyResult-method features features,ClassifyResult-method
#' models models,ClassifyResult-method performance
#' performance,ClassifyResult-method tunedParameters
#' tunedParameters,ClassifyResult-method totalPredictions
#' totalPredictions,ClassifyResult-method
#' @docType class
#' 
#' @section Constructor:
#' \preformatted{ClassifyResult(characteristics, originalNames, originalFeatures,
#'               rankedFeatures, chosenFeatures, models, tunedParameters, predictions, actualOutcomes, modellingParams = NULL, finalModel = NULL)}
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
#' \item{\code{rankedFeatures}}{All features, from most to least important. Character vector
#' or a data frame if data set has multiple kinds of measurements on the same set of samples.}
#' \item{\code{chosenFeatures}}{Features selected at each fold. Character
#' vector or a data frame if data set has multiple kinds of measurements on the same set of samples.}
#' \item{\code{models}}{All of the models fitted to the training data.}
#' \item{\code{tunedParameters}}{Names of tuning parameters and the value chosen of each parameter.}
#' \item{\code{predictions}}{A data frame containing sample IDs, predicted class or risk and information about the 
#' cross-validation iteration in which the prediction was made.}
#' \item{\code{actualOutcomes}}{The known class or survival data of each sample.}
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
#' \item{\code{allFeatureNames(result)}}{Returns a vector of features present in the data set.}}
#' \describe{
#' \item{\code{actualOutcomes(result)}}{Returns the known outcomes of each sample.}}
#' \describe{
#' \item{\code{models(result)}}{A \code{list} of the models fitted for each training.}}
#' \describe{
#' \item{\code{chosenFeatureNames(result)}}{A \code{list} of the features selected for each training.}}
#' \describe{
#' \item{\code{predictions(result)}}{Returns a \code{data.frame} which has columns with test sample,
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
#' @usage NULL
#' @export
setGeneric("ClassifyResult", function(characteristics, originalNames, originalFeatures, ...)
standardGeneric("ClassifyResult"))

#' @rdname ClassifyResult-class
#' @exportClass ClassifyResult
setClass("ClassifyResult", representation(
  characteristics = "DataFrame",
  originalNames = "character",
  originalFeatures = "characterOrDataFrame",
  rankedFeatures = "listOrNULL",
  chosenFeatures = "listOrNULL",
  actualOutcomes = "factorOrSurv",
  models = "list",
  tune = "listOrNULL",
  predictions = "data.frame",
  performance = "listOrNULL",
  modellingParams = "ModellingParamsOrNULL",
  finalModel = "listOrNULL")
)
#' @rdname ClassifyResult-class
#' @usage NULL
#' @export
setMethod("ClassifyResult", c("DataFrame", "character", "characterOrDataFrame"),
          function(characteristics, originalNames, originalFeatures,
                   rankedFeatures, chosenFeatures, models, tunedParameters, predictions, actualOutcomes, modellingParams = NULL, finalModel = NULL)
          {
            new("ClassifyResult", characteristics = characteristics,
                originalNames = originalNames, originalFeatures = originalFeatures,
                rankedFeatures = rankedFeatures, chosenFeatures = chosenFeatures,
                models = models, tune = tunedParameters,
                predictions = predictions, actualOutcomes = actualOutcomes, modellingParams = modellingParams, finalModel = finalModel)
          })

#' @usage NULL
#' @export
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
setMethod("sampleNames", c("ClassifyResult"),
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
setMethod("chosenFeatureNames", c("ClassifyResult"),
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
setMethod("models", c("ClassifyResult"),
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
setMethod("predictions", c("ClassifyResult"),
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
setMethod("performance", c("ClassifyResult"),
          function(object)
          {
            object@performance
          })

#' @export
#' @usage NULL
setGeneric("actualOutcomes", function(object, ...)
standardGeneric("actualOutcomes"))

#' @rdname ClassifyResult-class
#' @usage NULL
#' @export
setMethod("actualOutcomes", c("ClassifyResult"),
          function(object)
          {
            object@actualOutcomes
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
setMethod("totalPredictions", c("ClassifyResult"),
          function(result)
          {
              nrow(predictions(result))
          })

setClass("MixModelsListsSet", representation(set = "list"))
setGeneric("MixModelsListsSet", function(set, ...) standardGeneric("MixModelsListsSet"))
setMethod("MixModelsListsSet", c("list"), function(set) new("MixModelsListsSet", set = set))