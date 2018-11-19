setGeneric("pairsDifferencesSelection", function(measurements, ...)
           {standardGeneric("pairsDifferencesSelection")})

setMethod("pairsDifferencesSelection", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, featurePairs = NULL, ...)
{
  if(is.null(featurePairs))
    stop("No feature pairs provided but some must be.")
  if(!"Pairs" %in% class(featurePairs))
    stop("'featurePairs' must be of type Pairs.")
            
  pairsDifferencesSelection(DataFrame(t(measurements), check.names = FALSE), classes, featurePairs, ...)
})

setMethod("pairsDifferencesSelection", "DataFrame",
          function(measurements, classes, featurePairs = NULL, datasetName,
                   trainParams, predictParams, resubstituteParams,
                   selectionName = "Pairs Differences", verbose = 3)
{
  if(is.null(featurePairs))
    stop("No feature pairs provided but some must be.")
  if(!"Pairs" %in% class(featurePairs))
    stop("'featurePairs' must be of type Pairs.")
            
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  
  if(!all(S4Vectors::first(featurePairs) %in% colnames(measurements)) && !all(S4Vectors::second(featurePairs) %in% colnames(measurements)))
    stop("Some interactors are not found in 'measurements'. Ensure that 'featurePairs' does not have
         any features not in 'measurements'.")  
  
  if(verbose == 3)
    message("Selecting pairs of features with consistent differences.")


  oneClassTraining <- which(classes == levels(classes)[1])
  otherClassTraining <- which(classes == levels(classes)[2])
  oneClassMeasurements <- measurements[oneClassTraining, ]
  otherClassMeasurements <- measurements[otherClassTraining, ]
  
  oneClassDifferences <- sapply(1:length(featurePairs), function(pairIndex)
                         {
                           sum(oneClassMeasurements[, S4Vectors::first(featurePairs[pairIndex])] - oneClassMeasurements[, S4Vectors::second(featurePairs[pairIndex])])
                         })
  
  otherClassDifferences <- sapply(1:length(featurePairs), function(pairIndex)
                           {
                             sum(otherClassMeasurements[, S4Vectors::first(featurePairs[pairIndex])] - otherClassMeasurements[, S4Vectors::second(featurePairs[pairIndex])])
                           })
  pairsClassDifferences <- otherClassDifferences - oneClassDifferences
  orderedFeatures <- order(abs(pairsClassDifferences), decreasing = TRUE)

  .pickFeatures(measurements, classes, featurePairs, datasetName, trainParams, predictParams,
                resubstituteParams, orderedFeatures, selectionName, verbose)  
})

# One or more omics data sets, possibly with clinical data.
setMethod("pairsDifferencesSelection", "MultiAssayExperiment",
          function(measurements, target = names(measurements)[1], featurePairs = NULL, ...)
{
  if(is.null(featurePairs))
    stop("No feature pairs provided but some must be.")
  if(!"Pairs" %in% class(featurePairs))
    stop("'featurePairs' must be of type Pairs.")         
            
  tablesAndClasses <- .MAEtoWideTable(measurements, target)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]            
  pairsDifferencesSelection(dataTable, classes, featurePairs, ...)
})