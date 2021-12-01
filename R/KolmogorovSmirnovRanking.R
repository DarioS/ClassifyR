setGeneric("KolmogorovSmirnovRanking", function(measurements, ...)
           standardGeneric("KolmogorovSmirnovRanking"))

# Matrix of numeric measurements.
setMethod("KolmogorovSmirnovRanking", "matrix", function(measurements, classes, ...)
{
  KolmogorovSmirnovRanking(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("KolmogorovSmirnovRanking", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, ..., verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")

  if(verbose == 3)
    message("Selecting features by Kolmogorov Smirnov distance.")

  oneClass <- classes == levels(classes)[1]
  otherClass <- classes == levels(classes)[2]
  KSdistance <- apply(measurements, 2, function(featureColumn)
                      stats::ks.test(featureColumn[oneClass], featureColumn[otherClass], ...)[["statistic"]])

  order(KSdistance, decreasing = TRUE)
})

# One or more omics data sets, possibly with clinical data.
setMethod("KolmogorovSmirnovRanking", "MultiAssayExperiment",
function(measurements, targets = names(measurements), ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
            
  if(ncol(dataTable) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    KolmogorovSmirnovRanking(dataTable, classes, ...)
})