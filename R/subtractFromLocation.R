setGeneric("subtractFromLocation", function(measurements, ...)
           {standardGeneric("subtractFromLocation")})

setMethod("subtractFromLocation", "matrix", 
          function(measurements, training, location = c("mean", "median"),
                   absolute = TRUE, verbose = 3)
{
  location <- match.arg(location)
  measurementsTrain <- measurements[, training]
  if(location == "mean")
    featureTrainingLocations <- rowMeans(measurementsTrain, na.rm = TRUE)
  else # median.
    featureTrainingLocations <- apply(measurementsTrain, 1, median, na.rm = TRUE)
  transformed <- apply(measurements, 2, '-', featureTrainingLocations)
  if(absolute == TRUE)
    transformed <- abs(transformed)
  
  if(verbose == 3)
    message("Subtraction from ", location,
            {if(absolute == TRUE) " and absolute transformation"}, " completed.")
            
  transformed # Return an endomorphic variable; a matrix.
})

setMethod("subtractFromLocation", "DataFrame", 
          function(measurements, training, location = c("mean", "median"),
                   absolute = TRUE, verbose = 3)
          {
            whichNumeric <- sapply(measurements, is.numeric)
            if(sum(whichNumeric) == 0)
              stop("All features are not numeric but at least one must be.")
            
            if(sum(whichNumeric) != ncol(measurements) && verbose == 3)
              message("Some columns are not numeric. Only subtracting values from location for\n", 
                      "the columns containing numeric data.")
            
            location <- match.arg(location)
            measurementsTrain <- measurements[training, whichNumeric]
            if(location == "mean")
              featureTrainingLocations <- colMeans(measurementsTrain, na.rm = TRUE)
            else # median.
              featureTrainingLocations <- apply(measurementsTrain, 2, median, na.rm = TRUE)
            transformed <- measurements
            transformed[, whichNumeric] <- apply(measurements[, whichNumeric], 1, '-', featureTrainingLocations)
            if(absolute == TRUE)
              transformed <- abs(transformed)
            
            if(verbose == 3)
              message("Subtraction from ", location,
                      {if(absolute == TRUE) " and absolute transformation"}, " completed.")
            
            transformed # Return an endomorphic variable; a DataFrame.
          })

setMethod("subtractFromLocation", "MultiAssayExperiment", 
          function(measurements, training, targets = names(measurements),
                   location = c("mean", "median"), absolute = TRUE, verbose = 3)
{
  location <- match.arg(location)
  if(!all(targets %in% c(names(measurements), "colData")))
    stop("Some table names in 'targets' are not assay names in 'measurements' or \"colData\".")  
  
  if("colData" %in% targets)
  {
    targets <- targets[-which(match, "colData")]
    whichNumeric <- sapply(colData(measurements), is.numeric)
    clinicalTrain <- colData(measurements)[training, whichNumeric]
    if(location == "mean")
      featureTrainingLocations <- colMeans(clinicalTrain, na.rm = TRUE)
    else # median.
      featureTrainingLocations <- apply(clinicalTrain, 2, median, na.rm = TRUE)
    transformedClinical <- apply(colData(measurements), 1, '-', locations)
    if(absolute == TRUE)
      transformedClinical <- abs(transformedClinical)
    colData(measurements)[, whichNumeric] <- transformedClinical
  }
  
  measurementsTrain <- measurements[, training, targets]
  if(location == "mean")
    featureTrainingLocations <- lapply(experiments(measurementsTrain), rowMeans, na.rm = TRUE)
  else # median.
    featureTrainingLocations <- lapply(experiments(measurementsTrain), function(dataTable)
                                                              apply(dataTable, 1, median, na.rm = TRUE))
  transformedTables <- mapply(function(dataTable, locations) apply(dataTable, 2, '-', locations),
                              experiments(measurements[, , targets]), featureTrainingLocations, SIMPLIFY = FALSE)
  if(absolute == TRUE)
    transformedTables <- lapply(transformedTables, abs)
  
  transformed <- measurements
  experiments(transformed)[targets] <- ExperimentList(transformedTables)

  if(verbose == 3)
    message("Subtraction from ", location,
            {if(absolute == TRUE) " and absolute transformation"}, " completed.")
  
  transformed # Return an endomorphic variable; a MultiAssayExperiment.
})