#' Perform a Single Classification
#' 
#' For a data set of features and samples, the classification process is run.
#' It consists of data transformation, feature selection, classifier training
#' and testing.
#' 
#' This function only performs one classification and prediction. See
#' \code{\link{runTests}} for a driver function that enables a number of
#' different cross-validation schemes to be applied and uses this function to
#' perform each iteration.
#' 
#' @aliases runTest runTest,matrix-method runTest,DataFrame-method
#' runTest,MultiAssayExperiment-method
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
#' @param outcomesTrain Either a factor vector of classes, a \code{\link{Surv}} object, or
#' a character string, or vector of such strings, containing column name(s) of column(s)
#' containing either classes or time and event information about survival.
#' @param crossValParams An object of class \code{\link{CrossValParams}},
#' specifying the kind of cross-validation to be done, if nested
#' cross-validation is used to tune any parameters.
#' @param modellingParams An object of class \code{\link{ModellingParams}},
#' specifying the class rebalancing, transformation (if any), feature selection
#' (if any), training and prediction to be done on the data set.
#' @param targets If \code{measurementsTrain} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"sampleInfo"} is also a valid value
#' and specifies that numeric variables from the sample information data table will be
#' used.
#' @param outcomeColumn The column name of the sample information data table containing
#' the outcome to be predicted. Will automatically be removed from the sample information
#' data table during model training. Must be specified if data is a
#' \code{MultiAssayExperiment}.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method.
#' @param characteristics A \code{\link{DataFrame}} describing the
#' characteristics of the classification used. First column must be named
#' \code{"charateristic"} and second column must be named \code{"value"}.
#' Useful for automated plot annotation by plotting functions within this
#' package.  Transformation, selection and prediction functions provided by
#' this package will cause the characteristics to be automatically determined
#' and this can be left blank.
#' @param training A vector which specifies the training samples.
#' @param testing A vector which specifies the test samples.
#' @param verbose Default: 1. A number between 0 and 3 for the amount of
#' progress messages to give.  A higher number will produce more messages as
#' more lower-level functions print messages.
#' @param .iteration Not to be set by a user. This value is used to keep track
#' of the cross-validation iteration, if called by \code{\link{runTests}}.
#' @return If called directly by the user rather than being used internally by
#' \code{\link{runTests}}, a \code{\link{ClassifyResult}} object. Otherwise a
#' list of different aspects of the result which is passed back to
#' \code{\link{runTests}}.
#' @author Dario Strbenac
#' @examples
#' 
#'   #if(require(sparsediscrim))
#'   #{
#'     data(asthma)
#'     tuneList <- list(nFeatures = seq(5, 25, 5), performanceType = "Balanced Error")
#'     selectParams <- SelectParams(limmaRanking, tuneParams = tuneList)
#'     modellingParams <- ModellingParams(selectParams = selectParams)
#'     runTest(measurements, classes, modellingParams = modellingParams,
#'             training = seq(1, ncol(measurements), 2), testing = seq(2, ncol(measurements), 2))
#'   #}
#' 
#' @importFrom S4Vectors do.call
#' @export
setGeneric("runTest", function(measurementsTrain, ...)
           standardGeneric("runTest"))

setMethod("runTest", "matrix", # Matrix of numeric measurements.
  function(measurementsTrain, outcomesTrain, ...)
{
  if(is.null(rownames(measurementsTrain)))
    stop("'measurementsTrain' matrix must have sample identifiers as its row names.")    
  runTest(DataFrame(measurementsTrain, check.names = FALSE), outcomesTrain, ...)
})

setMethod("runTest", "DataFrame", # Sample information data or one of the other inputs, transformed.
function(measurementsTrain, outcomesTrain, measurementsTest, outcomesTest,
         crossValParams = CrossValParams(), # crossValParams used for tuning optimisation.
         modellingParams = ModellingParams(), characteristics = DataFrame(), verbose = 1, .iteration = NULL)
{
  if(is.null(.iteration)) # Not being called by runTests but by user. So, check the user input.
  {
    if(is.null(rownames(measurementsTrain)))
      stop("'measurementsTrain' DataFrame must have sample identifiers as its row names.")
    if(any(is.na(measurementsTrain)))
      stop("Some data elements are missing and classifiers don't work with missing data. Consider imputation or filtering.")                
    
    splitDatasetTrain <- .splitDataAndOutcomes(measurementsTrain, outcomesTrain)
    # Rebalance the class sizes of the training samples by either downsampling or upsampling
    # or leave untouched if balancing is none.
    if(!is(classes, "Surv"))
    {
      rebalancedTrain <- .rebalanceTrainingClasses(splitDatasetTrain[["measurements"]], splitDatasetTrain[["outcomes"]], training, testing, modellingParams@balancing)
      measurementsTrain <- rebalanced[["measurementsTrain"]]
      classesTrain <- rebalanced[["classesTrain"]]
    }
  }
  
  # All input features.
  if(!is.null(S4Vectors::mcols(measurementsTrain)))
    allFeatures <- S4Vectors::mcols(measurementsTrain)
  else
    allFeatures <- colnames(measurementsTrain)
  
  if(!is.null(modellingParams@transformParams))
  {
    if(length(modellingParams@transformParams@intermediate) != 0)
      modellingParams@transformParams <- .addIntermediates(modellingParams@transformParams)
    
    measurementsTrain <- tryCatch(.doTransform(measurementsTrain, modellingParams@transformParams, verbose), error = function(error) error[["message"]])
    if(is.character(measurementsTrain)) return(measurementsTrain) # An error occurred.
  }

  rankedFeatures <- NULL
  selectedFeatures <- NULL
  tuneDetailsSelect <- NULL
  if(!is.null(modellingParams@selectParams))
  {
    if(length(modellingParams@selectParams@intermediate) != 0)
      modellingParams@selectParams <- .addIntermediates(modellingParams@selectParams)
    
    topFeatures <- tryCatch(.doSelection(measurementsTrain, outcomesTrain, modellingParams, verbose = verbose, crossValParams = crossValParams),
                            error = function(error) error[["message"]]) 

    if(is.character(topFeatures)) return(topFeatures) # An error occurred.
    rankedFeatures <- topFeatures[[1]] # Extract for result object.
    selectedFeatures <- topFeatures[[2]] # Extract for subsetting.
    tuneDetailsSelect <- topFeatures[[3]]
  
    if(modellingParams@selectParams@subsetToSelections == TRUE)
    { # Subset the the data table to only the selected features.
      if(is.null(S4Vectors::mcols(measurementsTrain)))
      { # Input was ordinary matrix or DataFrame.
          measurementsTrain <- measurementsTrain[, selectedFeatures, drop = FALSE]
      } else { # Input was MultiAssayExperiment. # Match the selected features to the data frame columns
          selectedIDs <-  do.call(paste, selectedFeatures)
          featuresIDs <- do.call(paste, S4Vectors::mcols(measurementsTrain)[, c("dataset", "feature")])
          selectedColumns <- match(selectedIDs, featuresIDs)
          measurementsTrain <- measurementsTrain[, selectedColumns, drop = FALSE]
      }
    }
  } 
  
  # Training stage.
  if(length(modellingParams@trainParams@intermediate) > 0)
    modellingParams@trainParams <- .addIntermediates(modellingParams@trainParams)
  
  trained <- tryCatch(.doTrain(measurementsTrain, classes, training, testing, modellingParams, verbose),
                      error = function(error) error[["message"]])

  if(is.character(trained)) return(trained) # An error occurred.
    
  tuneDetailsTrain <- trained[[2]] # Second element is tuning results.
  
  
  if(!is.null(modellingParams@trainParams@getFeatures)) # Features chosen inside classifier.
  {
    extrasList <- list()
    extras <- .methodFormals(modellingParams@trainParams@getFeatures)[-1]
    if(length(extras) > 0)
      extrasList <- mget(names(extras))
    
    featureInfo <- do.call(modellingParams@trainParams@getFeatures, c(trained[[1]], extrasList))
    rankedFeatures <- featureInfo[[1]]
    selectedFeatures <- featureInfo[[2]]
  }
  
  if(!is.null(modellingParams@predictParams))
  {
    if(length(modellingParams@predictParams@intermediate) != 0)
      modellingParams@predictParams <- .addIntermediates(modellingParams@predictParams)
                             
    predictedClasses <- tryCatch(.doTest(trained[["model"]], measurementsTrain, testing, modellingParams@predictParams, verbose),
                                error = function(error) error[["message"]]
                                )
    if(is.character(predictedClasses)) # An error occurred.
      return(predictedClasses) # Return early.
    
  } else {
    predictedClasses <- trained[[1]]
  }
  
  if(is.null(modellingParams@predictParams)) models <- NULL else models <- trained[[1]] # One function for training and testing. Typically, the models aren't returned to the user, such as Poisson LDA implemented by PoiClaClu.
  if(!is.null(tuneDetailsSelect)) tuneDetails <- tuneDetailsSelect else tuneDetails <- tuneDetailsTrain
  if(!is.null(.iteration)) # This function was called by runTests.
  {
    list(ranked = rankedFeatures, selected = selectedFeatures, models = models, testSet = testingSamplesIDs, predictions = predictedClasses, tune = tuneDetails)
  } else { # runTest is being used directly, rather than from runTests. Create a ClassifyResult object.
    # Only one training, so only one tuning choice, which can be summarised in characteristics.
    modParamsList <- list(modellingParams@transformParams, modellingParams@selectParams, modellingParams@trainParams, modellingParams@predictParams)
    if(!is.null(tuneDetails)) characteristics <- rbind(characteristics, data.frame(characteristic = colnames(tuneDetails),
                                                                                   value = unlist(tuneDetails)))
    autoCharacteristics <- do.call(rbind, lapply(modParamsList, function(stageParams) if(!is.null(stageParams)) stageParams@characteristics))
    characteristics <- .filterCharacteristics(characteristics, autoCharacteristics)
    characteristics <- rbind(characteristics, S4Vectors::DataFrame(characteristic = "Cross-validation", value = "Independent Set"))
.
    extras <- lapply(modParamsList, function(stageParams) if(!is.null(stageParams))stageParams@otherParams)
    extrasDF <- DataFrame(characteristic = names(extras), value = unlist(extras))
    characteristics <- rbind(characteristics, extrasDF)
    
    ClassifyResult(characteristics, rownames(measurementsTrain), allFeatures, list(rankedFeatures), list(selectedFeatures),
                   list(models), tuneDetails, data.frame(sample = testing, class = predictedClasses), classes)
  }  
})

setMethod("runTest", c("MultiAssayExperiment"),
          function(measurementsTrain, measurementsTest, targets = names(measurements), outcomesColumns, ...)
{
  tablesAndClassesTrain <- .MAEtoWideTable(measurementsTrain, targets, outcomesColumns, restrict = NULL)
  tablesAndClassesTest <- .MAEtoWideTable(measurementsTest, targets, outcomesColumns, restrict = NULL)
  runTest(tablesAndClassesTrain[["dataTable"]], tablesAndClassesTrain[["outcomes"]],
          tablesAndClassesTest[["dataTable"]], tablesAndClassesTest[["outcomes"]], ...)            
})