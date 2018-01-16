setGeneric("bartlettSelection", function(measurements, ...)
{standardGeneric("bartlettSelection")})

setMethod("bartlettSelection", "matrix",
          function(measurements, classes, ...)
          {
            .bartlettSelection(DataFrame(t(measurements), check.names = FALSE), classes, ...)
          })

setMethod("bartlettSelection", "DataFrame", # Clinical data only.
          function(measurements, classes, ...)
          {
            splitDataset <- .splitDataAndClasses(measurements, classes)
            measurements <- splitDataset[["measurements"]]
            isNumeric <- apply(measurements, 2, is.numeric)
            measurements <- measurements[, isNumeric, drop = FALSE]
            if(sum(isNumeric) == 0)
              stop("All features are not numeric but at least one must be.")
            .bartlettSelection(measurements, splitDataset[["classes"]], ...)
          })

# One or more omics datasets, possibly with clinical data.
setMethod("bartlettSelection", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), ...)
          {
            if(!all(targets %in% c(names(measurements), "colData")))
              stop("Some table names in 'targets' are not assay names in 'measurements' or \"colData\".")  
            
            measurements <- measurements[, , targets]
            allFeaturesMatrix <- wideFormat(measurements, colDataCols = colnames(colData(measurements)),
                                            check.names = FALSE)
            classes <- allFeaturesMatrix[, "class"]
            allFeaturesMatrix <- allFeaturesMatrix[, -match(c("primary", "class"),
                                                           colnames(allFeaturesMatrix))]
            isNumeric <- apply(allFeaturesMatrix, 2, is.numeric)
            allFeaturesMatrix <- allFeaturesMatrix[, isNumeric, drop = FALSE]
            if(ncol(allFeaturesMatrix) == 0)
              stop("No variables in data tables specified by \'targets\' are numeric.")
            else
              .bartlettSelection(allFeaturesMatrix, classes, ...)
          })

.bartlettSelection <- function(measurements, classes,
                               datasetName, trainParams, predictParams, resubstituteParams,
                               selectionName = "Bartlett Test", verbose = 3)
{
  if(verbose == 3)
    message("Calculating Bartlett statistic for each feature.")

  pValues <- apply(measurements, 2, function(featureColumn)
    bartlett.test(featureColumn, classes)[["p.value"]])
  orderedFeatures <- order(pValues)
  
  .pickFeatures(measurements, classes, datasetName, trainParams, predictParams,
                resubstituteParams, orderedFeatures, selectionName, verbose)
}