setGeneric("featureSetSummary", function(measurements, ...)
           {standardGeneric("featureSetSummary")})

setMethod("featureSetSummary", "matrix", # Matrix of numeric measurements.
          function(measurements, location = c("median", "mean"),
                   featureSets, minimumOverlapPercent = 80, verbose = 3)
{
  if(class(featureSets) != "FeatureSetCollection")
    stop("'featureSets' is not of type FeatureSetCollection but must be.")

  assayedFeatures <- rownames(measurements)
  keepSets <- sapply(featureSets@sets, function(featureSet)
    length(intersect(featureSet, assayedFeatures)) / length(featureSet) * 100 > minimumOverlapPercent)
  if(all(keepSets == FALSE))
    stop("No feature sets had an overlap of at least ", minimumOverlapPercent,
         "% with the data set's feature identifiers.")
  
  if(any(keepSets == FALSE)) # Filter out those sets without adequate identifier overlap.
  {
    if(verbose == 3)
      message("Based on", paste(minimumOverlapPercent, "% overlap rule, reducing", sep = ''), length(featureSets), "feature sets to", sum(keepSets), "feature sets.")
    featureSets <- featureSets[keepSets]
  }
  
  # Reduce set representations to only those features which were assayed.
  featureSets <- lapply(featureSets@sets, function(featureSet) intersect(featureSet, assayedFeatures))
  
  location <- match.arg(location)
  if(location == "mean")
    locationFunction <- mean
  else
    locationFunction <- median
  
    if(verbose == 3)
      message("Summarising features to feature sets.")

  # Transform measurements into one feature per set.
  apply(measurements, 2, function(sampleMeasurements)
  {
    sapply(featureSets, function(featureSet) locationFunction(sampleMeasurements[featureSet]))
  })
})

setMethod("featureSetSummary", "DataFrame", # Possibly mixed data types.
          function(measurements, location = c("median", "mean"),
                   featureSets, minimumOverlapPercent = 80, verbose = 3)
{
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  if(class(featureSets) != "FeatureSetCollection")
    stop("'featureSets' is not of type FeatureSetCollection but must be.")

  assayedFeatures <- colnames(measurements)
  keepSets <- sapply(featureSets@sets, function(featureSet)
    length(intersect(featureSet, assayedFeatures)) / length(featureSet) * 100 > minimumOverlapPercent)
  if(all(keepSets == FALSE))
    stop("No feature sets had an overlap of at least ", minimumOverlapPercent,
         "% with the data set's feature identifiers.")
  
  if(any(keepSets == FALSE)) # Filter out those sets without adequate identifier overlap.
  {
    if(verbose == 3)
      message("Based on", paste(minimumOverlapPercent, "% overlap rule, reducing", sep = ''), length(featureSets), "feature sets to", sum(keepSets), "feature sets.")
    featureSets <- featureSets[keepSets]
  }
  
  # Reduce set representations to only those features which were assayed.
  featureSets <- lapply(featureSets@sets, function(featureSet) intersect(featureSet, assayedFeatures))
  
  location <- match.arg(location)
  if(location == "mean")
    locationFunction <- mean
  else
    locationFunction <- median
  
    if(verbose == 3)
      message("Summarising features to feature sets.")

  # Transform measurements into one feature per set.
  measurements <- as.matrix(measurements)
  measurementsCollapsed <- t(apply(measurements, 1, function(sampleMeasurements)
  {
    sapply(featureSets, function(featureSet) locationFunction(sampleMeasurements[featureSet]))
  }))

  DataFrame(measurementsCollapsed, check.names = FALSE)
})

setMethod("featureSetSummary", "MultiAssayExperiment", # Pick one numeric table from the data set.
          function(measurements, target = NULL, location = c("median", "mean"),
                   featureSets, minimumOverlapPercent = 80, verbose = 3)
{
  if(is.null(target))
    stop("'target' is NULL but must specify one of the data sets in 'measurements'.")
  if(class(featureSets) != "FeatureSetCollection")
    stop("'featureSets' is not of type FeatureSetCollection but must be.")

  datasetUsed <- measurements[[target]]                                   
  assayedFeatures <- rownames(datasetUsed)
  keepSets <- sapply(featureSets@sets, function(featureSet)
    length(intersect(featureSet, assayedFeatures)) / length(featureSet) * 100 > minimumOverlapPercent)
  if(all(keepSets == FALSE))
    stop("No feature sets had an overlap of at least ", minimumOverlapPercent,
         "% with the data set's feature identifiers.")
  
  if(any(keepSets == FALSE)) # Filter out those sets without adequate identifier overlap.
  {
    if(verbose == 3)
      message("Based on", paste(minimumOverlapPercent, "% overlap rule, reducing", sep = ''), length(featureSets), "feature sets to", sum(keepSets), "feature sets.")
    featureSets <- featureSets[keepSets]
  }
  
  # Reduce set representations to only those features which were assayed.
  featureSets <- lapply(featureSets@sets, function(featureSet) intersect(featureSet, assayedFeatures))
  
  location <- match.arg(location)
  if(location == "mean")
    locationFunction <- mean
  else
    locationFunction <- median
  
    if(verbose == 3)
      message("Summarising features to feature sets.")

  # Transform measurements into one feature per set.
  transformed <- apply(datasetUsed, 2, function(sampleMeasurements)
  {
    sapply(featureSets, function(featureSet) locationFunction(sampleMeasurements[featureSet]))
  })
  
  measurements[[target]] <- transformed
  measurements
})