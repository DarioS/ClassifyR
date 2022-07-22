setGeneric("selectMulti", function(measurementsTrain, classesTrain, params, ...)
    standardGeneric("selectMulti"))

# DataFrame of numeric measurements, likely created by runTests or runTest.
setMethod("selectMulti", "DataFrame",
          function(measurementsTrain, classesTrain, params, verbose = 0)
          {
              assayTrain <- sapply(unique(mcols(measurementsTrain)[["Renamed Assay"]]), function(assay) measurementsTrain[, mcols(measurementsTrain)[["Renamed Assay"]] %in% assay], simplify = FALSE)
              featuresIndices <- mapply(.doSelection, 
                                         measurements = assayTrain,
                                         modellingParams = params,
                                         MoreArgs = list(outcomeTrain = classesTrain, 
                                                         crossValParams = CrossValParams(permutations = 1, folds = 5), ###### Where to get this from?
                                                         verbose = 0), SIMPLIFY = FALSE
                                        )
              
              unique(unlist(lapply(featuresIndices, "[[", 2)))
          })
