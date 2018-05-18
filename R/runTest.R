setGeneric("runTest", function(measurements, ...)
           {standardGeneric("runTest")})

setMethod("runTest", c("matrix"), # Matrix of numeric measurements.
  function(measurements, classes, ...)
{
  runTest(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("runTest", c("DataFrame"), # Clinical data only.
function(measurements, classes, datasetName, classificationName, training, testing,
         params = list(SelectParams(), TrainParams(), PredictParams()),
         verbose = 1, .iteration = NULL)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)

  stagesParamClasses <- sapply(params, class)
  if(match("TrainParams", stagesParamClasses) > match("PredictParams", stagesParamClasses))
    stop("\"PredictParams\" variable must not be before \"TrainParams\" in 'params' list.")

  transformParams <- params[[match("TransformParams", stagesParamClasses)]]
  selectParams <- params[[match("SelectParams", stagesParamClasses)]]
  trainParams <- params[[match("TrainParams", stagesParamClasses)]]
  predictParams <- params[[match("PredictParams", stagesParamClasses)]]
  
  if(!is.null(S4Vectors::mcols(measurements)))
    allFeatures <- S4Vectors::mcols(measurements)
  else
    allFeatures <- colnames(measurements)
  
  lastSize <- 1
  for(stageIndex in 1:length(params))
  {
    switch(stagesParamClasses[[stageIndex]],
           TransformParams = {
                               if(length(transformParams@intermediate) != 0)
                                 transformParams@otherParams <- c(transformParams@otherParams, mget(transformParams@intermediate))

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
                                                     if(is.null(S4Vectors::mcols(variety)) == TRUE)
                                                     { # Input was ordinary matrix or DataFrame.
                                                       variety[, selectedFeatures, drop = FALSE]
                                                     } else { # Input was MultiAssayExperiment.
                                                       selectedColumns <- apply(selectedFeatures, 2, function(selectedFeature)
                                                       {
                                                         intersect(which(selectedFeature[1] == S4Vectors::mcols(variety)[, "dataset"]),
                                                                   which(selectedFeature[2] == S4Vectors::mcols(variety)[, "feature"]))
                                                       })
                                                       variety <- variety[, selectedColumns, drop = FALSE]
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
                                                         if(is.null(S4Vectors::mcols(variety)) == TRUE)
                                                         { # Input was ordinary matrix or DataFrame.
                                                           variety[, features, drop = FALSE]
                                                         } else { # Input was MultiAssayExperiment.
                                                           selectedColumns <- apply(selectedFeatures, 2, function(selectedFeature)
                                                           {
                                                             intersect(which(selectedFeature[1] == S4Vectors::mcols(variety)[, "dataset"]),
                                                                       which(selectedFeature[2] == S4Vectors::mcols(variety)[, "feature"]))
                                                           })
                                                           variety <- variety[, selectedColumns, drop = FALSE]
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

                              newSize <- if("list" %in% class(trained)) length(trained) else 1
                              if(newSize / lastSize != 1) # More varieties were created.
                              {
                                measurements <- unlist(lapply(if(class(measurements) == "list") measurements else list(measurements), function(variety)
                                                                                      lapply(1:(newSize / lastSize), function(replicate) variety)),
                                                                               recursive = FALSE)
                                names(measurements) <- names(trained)
                              }
                              
                              lastSize <- newSize
                              if("list" %in% class(trained))
                              {
                                tuneDetails <- lapply(trained, attr, "tune")
                                if(!is.null(trainParams@getFeatures)) # Features chosen inside classifier.
                                {
                                  featureInfo <- lapply(trained, trainParams@getFeatures)
                                  rankedFeatures <- lapply(featureInfo, '[[', 1)
                                  selectedFeatures <- lapply(featureInfo, '[[', 2)
                                }
                              } else {
                                tuneDetails <- attr(trained, "tune")
                                if(!is.null(trainParams@getFeatures)) # Features chosen inside classifier.
                                {                                
                                  rankedFeatures <- trainParams@getFeatures(trained)[[1]]
                                  selectedFeatures <- trainParams@getFeatures(trained)[[2]]
                                }
                              }
                                if(is.null(tuneDetails)) tuneDetails <- list(tuneDetails)
                              },
              PredictParams = {
                                 if(length(predictParams@intermediate) != 0)
                                   predictParams@otherParams <- c(predictParams@otherParams, mget(predictParams@intermediate))
                                  predictedClasses <- tryCatch(.doTest(trained, measurements, testing, predictParams, verbose),
                                                               error = function(error) error[["message"]])
                                  if(is.character(predictedClasses)) # An error occurred.
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
      return(ClassifyResult(datasetName, classificationName, selectParams@selectionName, rownames(measurements), allFeatures,
                            list(rankedFeatures), list(selectedFeatures), list(data.frame(sample = testing, class = predictedClasses)),
                            classes, list("independent"), tuneDetails)
             )
    } else { # A variety of predictions were made.
      if(!"list" %in% class(selectedFeatures))
      {
        rankedFeatures <- list(rankedFeatures)
        selectedFeatures <- list(selectedFeatures)
      }
      return(mapply(function(varietyPredictions, varietyTunes, varietyRanked, varietySelected)
      {
        if(is.null(varietyTunes)) varietyTunes <- list(varietyTunes)
        ClassifyResult(datasetName, classificationName, selectParams@selectionName, rownames(measurements), allFeatures,
                       list(varietyRanked), list(varietySelected), list(data.frame(sample = testing, class = varietyPredictions)),
                       classes, list("independent"), varietyTunes)
      }, predictedClasses, tuneDetails, rankedFeatures, selectedFeatures, SIMPLIFY = FALSE))
    }
  }  
})

setMethod("runTest", c("MultiAssayExperiment"),
          function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, restrict = NULL)
  runTest(tablesAndClasses[["dataTable"]], tablesAndClasses[["classes"]], ...)            
})