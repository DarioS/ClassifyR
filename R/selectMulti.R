setGeneric("selectMulti", function(measurementsTrain, classesTrain, params, ...)
    standardGeneric("selectMulti"))

# DataFrame of numeric measurements, likely created by runTests or runTest.
setMethod("selectMulti", "DataFrame",
          function(measurementsTrain, classesTrain, params, verbose = 0)
          {
              assayTrain <- sapply(unique(mcols(measurementsTrain)[["Renamed Dataset"]]), function(x) measurementsTrain[, mcols(measurementsTrain)[["Renamed Dataset"]] %in% x], simplify = FALSE)

              featuresIndices <- mapply(.doSelection, 
                                         measurements = assayTrain,
                                         modellingParams = params,
                                         MoreArgs = list(outcomesTrain = classesTrain, 
                                                         crossValParams = CrossValParams(permutations = 1, folds = 5), ###### Where to get this from?
                                                         verbose = 0), SIMPLIFY = FALSE
                                        )
              
              unique(unlist(lapply(featuresIndices, "[[", 2)))
          })
