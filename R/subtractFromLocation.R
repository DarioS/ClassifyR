setGeneric("subtractFromLocation", function(measurements, ...)
           {standardGeneric("subtractFromLocation")})

setMethod("subtractFromLocation", "matrix", 
          function(measurements, ...)
{
  groupsTable <- data.frame(row.names = names(classes)) # Classes are irrelevant.     
  subtractFromLocation(MultiAssayExperiment(ExperimentList(list(dataTable = measurements)),
                                            S4Vectors::DataFrame(groupsTable)), ...)
})

setMethod("subtractFromLocation", "MultiAssayExperiment", 
          function(measurements, training, targets = names(measurements),
                   location = c("mean", "median"), absolute = TRUE, verbose = 3)
{
  location <- match.arg(location)
  measurementsTrain <- measurements[, training, targets]
  if(location == "mean")
    featureTrainingLocations <- lapply(experiments(measurementsTrain), rowMeans)
  else # median.
    featureTrainingLocations <- lapply(experiments(measurementsTrain), function(dataTable)
                                                              apply(dataTable, 1, median))
  transformed <- mapply(function(dataTable, locations) apply(dataTable, 2, '-', locations),
                        experiments(measurements[, , targets]), featureTrainingLocations, SIMPLIFY = FALSE)
  if(absolute == TRUE)
    transformed <- lapply(transformed, abs)
  experiments(measurements)[targets] <- ExperimentList(transformed)
  if(verbose == 3)
    message("Subtraction from ", location,
            {if(absolute == TRUE) " and absolute transformation"}, " completed.")
  
  measurements
})