setGeneric("runTest", function(measurements, ...)
           {standardGeneric("runTest")})

setMethod("runTest", c("matrix"), # Matrix of numeric measurements.
  function(measurements, classes, ...)
{
  .runTest(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("runTest", c("DataFrame"), # Clinical data only.
function(measurements, classes, ...)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  .runTest(measurements, splitDataset[["classes"]], ...)
})

setMethod("runTest", c("MultiAssayExperiment"),
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, restrict = NULL)
  .runTest(tablesAndClasses[["dataTable"]], tablesAndClasses[["classes"]], ...)            
})

.runTest <- function(measurements, classes, datasetName, classificationName, training, testing,
                     params = list(SelectParams(), TrainParams(), PredictParams()),
                     verbose = 1, .iteration = NULL)
{
  stagesParamClasses <- sapply(params, class)
  if(match("TrainParams", stagesParamClasses) > match("PredictParams", stagesParamClasses))
    stop("\"PredictParams\" variable must not be before \"TrainParams\" in 'params' list.")

  transformParams <- params[[match("TransformParams", stagesParamClasses)]]
  selectParams <- params[[match("SelectParams", stagesParamClasses)]]
  trainParams <- params[[match("TrainParams", stagesParamClasses)]]
  predictParams <- params[[match("PredictParams", stagesParamClasses)]]
  
  allSamples <- switch(class(measurements), matrix = colnames(measurements),
                                            DataFrame = rownames(measurements),
                                            MultiAssayExperiment = rownames(colData(measurements)))
  
  allFeatures <- switch(class(measurements), matrix = rownames(measurements),
                                             DataFrame = colnames(measurements),
                                             MultiAssayExperiment = mcols(colData(measurements)))
  
  
  lastSize <- 1
  for(stageIndex in 1:length(params))
  {
    switch(stagesParamClasses[[stageIndex]],
           TransformParams = {
                               if(length(transformParams@intermediate) != 0)
                                 transformParams@otherParams <- c(transformParams@otherParams, mget(transformParams@intermediate))

                               transformParams@otherParams <- c(transformParams@otherParams, list(training = training))
                               measurements <- tryCatch(.doTransform(measurements, transformParams, verbose), error = function(error) error[["message"]])
                               if(is.character(measurements)) return(measurements) # An error occurred.
                               newSize <- if(class(measurements) == "list") length(measurements) else 1
                             },
              SelectParams = {
                               if(length(selectParams@intermediate) != 0)
                                 selectParams@otherParams <- c(selectParams@otherParams, mget(selectParams@intermediate))

                               topFeatures <- tryCatch(.doSelection(measurements, classes, training, selectParams,
                                                                trainParams, predictParams, verbose), error = function(error) error[["message"]])
                               if(is.character(topFeatures)) return(topFeatures) # An error occurred.

                               if(class(topFeatures[[2]]) == "list") # Check the chosen features list element, because a ranking is not present for ensemble selection.
                               {
                                 multiSelection <- TRUE
                               } else {
                                 multiSelection <- FALSE
                               }

                               rankedFeatures <- topFeatures[[1]] # Extract for result object.
                               selectedFeatures <- topFeatures[[2]] # Extract for subsetting.

                               if(selectParams@subsetToSelections == TRUE)
                               {
                                 if(multiSelection == FALSE)
                                 {
                                   if(class(measurements) != "list") # Put into list.
                                     measurements <- list(measurements)
                                   measurements <- lapply(measurements, function(variety)
                                                   {
                                                     if(class(variety) == "matrix")
                                                     {
                                                       variety[selectedFeatures, ]
                                                     } else if(class(variety) == "DataFrame") {
                                                       variety[, selectedFeatures]
                                                     } else {
                                                       assaysFeatures <- subset(selectedFeatures, dataset != "clinical")
                                                       sampleInfoFeatures <- subset(selectedFeatures, dataset == "clinical")
                                                       variety <- variety[assaysFeatures[, "feature"], , assaysFeatures[, "dataset"]]
                                                       colData(variety) <- colData(variety)[sampleInfoFeatures[, "variable"]]
                                                       variety
                                                     }
                                                   })
                                   if(length(measurements) == 1 && class(measurements) == "list")  # Restore to original container type.
                                     measurements <- measurements[[1]]
                                   measurements
                                 } else {
                                   if(class(measurements) != "list") # Put into list.
                                     measurements <- list(measurements)
                                   
                                   measurements <- lapply(measurements, function(variety)
                                                   {
                                                     lapply(selectedFeatures, function(features)
                                                     {
                                                         if(class(variety) == "matrix")
                                                         {   
                                                           variety[features, ]
                                                         } else if(class(variety) == "DataFrame") {
                                                           variety[, features]
                                                         } else { # MultiAssayExperiment
                                                           assaysFeatures <- subset(features, dataset != "clinical")
                                                           sampleInfoFeatures <- subset(features, dataset == "clinical")
                                                           variety <- variety[assaysFeatures[, "feature"], , assaysFeatures[, "dataset"]]
                                                           colData(variety) <- colData(variety)[sampleInfoFeatures[, "variable"]]
                                                           variety
                                                         }
                                                       })
                                                     })
                                   if(length(measurements) == 1 && class(measurements) == "list")  # Restore to original container type.
                                     measurements <- measurements[[1]]
                                 }
                               } else {
                                 if(multiSelection == TRUE) # Multiple selection varieties. Replicate the experimental data.
                                 {
                                   if(class(measurements) != "list")
                                     measurements <- lapply(selectedFeatures, function(features) measurements)
                                   else
                                     measurements <- lapply(measurements, function(variety)
                                                          lapply(selectedFeatures, function(features) variety))
                                 }
                               }

                               if(class(measurements) == "list" && class(measurements[[1]]) == "list")
                               {
                                 oldNames <- sapply(measurements, names)
                                 newNames <- unlist(lapply(measurements, names))
                                 measurements <- unlist(measurements, recursive = FALSE)
                                 names(measurements) <- paste(rep(oldNames, each = length(measurements[[1]])), newNames, sep = ',')
                               }
                               
                               if(class(measurements) == "list") newSize <- length(measurements) else newSize <- 1
                               lastSize <- newSize
                             }, 
              TrainParams = {
                              if(length(trainParams@intermediate) != 0)
                                trainParams@otherParams <- c(trainParams@otherParams, mget(trainParams@intermediate))

                              trained <- tryCatch(.doTrain(measurements, classes, training, testing, trainParams, predictParams, verbose),
                                                  error = function(error) error[["message"]])
                              if(is.character(trained)) return(trained) # An error occurred.

                              newSize <- if(class(trained) == "list") length(trained) else 1
                              if(newSize / lastSize != 1) # More varieties were created.
                              {
                                measurements <- unlist(lapply(if(class(measurements) == "list") measurements else list(measurements), function(variety)
                                                                                      lapply(1:(newSize / lastSize), function(replicate) variety)),
                                                                               recursive = FALSE)
                                names(measurements) <- names(trained)
                              }
                              
                              lastSize <- newSize
                              if(class(trained) == "list")
                                tuneDetails <- lapply(trained, attr, "tune")
                              else
                                tuneDetails <- attr(trained, "tune")
                                if(is.null(tuneDetails)) tuneDetails <- list(tuneDetails)
                              },
              PredictParams = {
                                 if(length(predictParams@intermediate) != 0)
                                   predictParams@otherParams <- c(predictParams@otherParams, mget(predictParams@intermediate))
                                  predictedClasses <- tryCatch(.doTest(trained, measurements, testing, predictParams, verbose),
                                                               error = function(error) error[["message"]])
                                  if(is.character(predictedClasses) && grepl("^Error", predictedClasses)) # An error occurred.
                                    return(predictedClasses) # Return early. Don't make a ClassifyResult below.
                               }
           )
    
  }
  if(class(testing) == "logical") testing <- which(testing)
  # Rankings and selections might not be explicitly returned, such as for random forest classifier.
  if(!exists("rankedFeatures")) rankedFeatures <- NULL
  if(!exists("selectedFeatures")) selectedFeatures <- NULL
  
  if(!is.null(.iteration)) # This function was called by runTests.
  {
    list(ranked = rankedFeatures, selected = selectedFeatures, testSet = testing, predictions = predictedClasses, tune = tuneDetails)
  } else { # runTest is being used directly, rather than from runTests. Create a ClassifyResult object.
    
    if(class(predictedClasses) != "list")
    {
      return(ClassifyResult(datasetName, classificationName, selectParams@selectionName, allSamples, allFeatures,
                            list(rankedFeatures), list(selectedFeatures), list(data.frame(sample = testing, label = predictedClasses)),
                            classes, list("independent"), tuneDetails)
             )
    } else { # A variety of predictions were made.
      return(mapply(function(varietyPredictions, varietyTunes)
      {
        if(is.null(varietyTunes)) varietyTunes <- list(varietyTunes)
        ClassifyResult(datasetName, classificationName, selectParams@selectionName, allSamples, allFeatures,
                       list(rankedFeatures), list(selectedFeatures), list(data.frame(sample = testing, label = varietyPredictions)),
                       classes, list("independent"), varietyTunes)
      }, predictedClasses, tuneDetails, SIMPLIFY = FALSE))
    }
  }
}