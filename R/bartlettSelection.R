setGeneric("bartlettSelection", function(measurements, ...)
{standardGeneric("bartlettSelection")})

setMethod("bartlettSelection", "matrix", 
          function(measurements, classes, ...)
          {
            groupsTable <- data.frame(class = classes, row.names = names(classes))
            measurementsSet <- MultiAssayExperiment(list(dataTable = measurements),
                                                    S4Vectors::DataFrame(groupsTable)
                                                   )
            bartlettSelection(measurementsSet, ...)
          })

setMethod("bartlettSelection", "DataFrame", # Clinical data only.
          function(measurements, classes, ...)
          {

          })

setMethod("bartlettSelection", "MultiAssayExperiment", # One or more omics datasets, possibly with clinical data.
          function(measurements, targets = names(measurements),
                   datasetName, trainParams, predictParams, resubstituteParams,
                   selectionName = "Bartlett Test", verbose = 3)
          {
            measurements <- measurements[, , targets]
            
            if(!all(sapply(experiments(measurements), is.numeric)))
               stop("At least one target table is not numeric but each target table must be.")
            if(verbose == 3)
              message("Calculating Bartlett statistic for each feature.")
            browser()
            allFeaturesMatrix <- wideFormat(measurements, colDataCols = "class", check.names = FALSE)
            classes <- allFeaturesMatrix[, "class"]
            allFeaturesMatrix <- allFeaturesMatrix[, -match(c("primary", "class"),
                                                           colnames(allFeaturesMatrix))]
            pValues <- apply(allFeaturesMatrix, 2, function(featureRow)
                                        bartlett.test(featureRow, classes)[["p.value"]])
            orderedFeatures <- order(pValues)
            .pickFeatures(allFeaturesMatrix, datasetName, trainParams, predictParams,
                          resubstituteParams, orderedFeatures, selectionName, verbose)
          })