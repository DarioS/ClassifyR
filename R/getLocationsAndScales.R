setGeneric("getLocationsAndScales", function(measurements, ...)
           {standardGeneric("getLocationsAndScales")})

setMethod("getLocationsAndScales", "matrix", # Matrix of numeric measurements.
          function(measurements, ...)
{
  getLocationsAndScales(DataFrame(t(measurements), check.names = FALSE), ...)
})

setMethod("getLocationsAndScales", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, location = c("mean", "median"), scale = c("SD", "MAD", "Qn"))
{
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  location <- match.arg(location) 
  scale <- match.arg(scale)
   
  if(scale == "Qn" && !requireNamespace("robustbase", quietly = TRUE))
    stop("The package 'robustbase' could not be found. Please install it.")
 
  setNames(list(switch(location,
                       mean = apply(measurements, 2, mean),
                       median = apply(measurements, 2, median)),
                switch(scale,
                       SD = apply(measurements, 2, sd),
                       MAD = apply(measurements, 2, mad),
                       Qn = apply(measurements, 2, robustbase::Qn))),
                c(location, scale))
})

# One or more omics data sets, possibly with clinical data.
setMethod("getLocationsAndScales", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
{
  if(!all(targets %in% c(names(measurements), "clinical")))
    stop("Some table names in 'targets' are not assay names in 'measurements' or \"clinical\".")  
            
  combinedData <- .MAEtoWideTable(measurements, targets)
  if(class(combinedData) == "list")
    combinedData <- combinedData[["dataTable"]]
  getLocationsAndScales(combinedData, ...)            
})