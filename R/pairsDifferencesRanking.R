setGeneric("pairsDifferencesRanking", function(measurements, ...)
           standardGeneric("pairsDifferencesRanking"))

setMethod("pairsDifferencesRanking", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, featurePairs = NULL, ...)
{
  if(is.null(featurePairs))
    stop("No feature pairs provided but some must be.")
  if(!"Pairs" %in% class(featurePairs))
    stop("'featurePairs' must be of type Pairs.")
            
  pairsDifferencesRanking(DataFrame(t(measurements), check.names = FALSE), classes, featurePairs, ...)
})

setMethod("pairsDifferencesRanking", "DataFrame",
          function(measurements, classes, featurePairs = NULL, verbose = 3)
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
  
  suppliedPairs <- length(featurePairs)
  keepPairs <- S4Vectors::first(featurePairs) %in% colnames(measurements) & S4Vectors::second(featurePairs) %in% colnames(measurements)
  featurePairs <- featurePairs[keepPairs]
  if(verbose == 3)
    message(suppliedPairs, " pairs input and ", length(featurePairs), " pairs remain after filtering based on data set row names.")
  
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
  order(abs(pairsClassDifferences), decreasing = TRUE)
})

# One or more omics data sets, possibly with clinical data.
setMethod("pairsDifferencesRanking", "MultiAssayExperiment",
          function(measurements, target = names(measurements)[1], featurePairs = NULL, ...)
{
  if(is.null(featurePairs))
    stop("No feature pairs provided but some must be.")
  if(!"Pairs" %in% class(featurePairs))
    stop("'featurePairs' must be of type Pairs.")         
            
  tablesAndClasses <- .MAEtoWideTable(measurements, target)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]            
  pairsDifferencesRanking(dataTable, classes, featurePairs, ...)
})