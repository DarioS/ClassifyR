setGeneric("selectMulti", function(measurementsTrain, classesTrain, params, ...)
    standardGeneric("selectMulti"))

# DataFrame of numeric measurements, likely created by runTests or runTest.
setMethod("selectMulti", "DataFrame",
          function(measurementsTrain, classesTrain, params, verbose = 0)
          {
              
              assayTrain <- sapply(unique(mcols(measurementsTrain)[["dataset"]]), function(x) measurementsTrain[,mcols(measurementsTrain)[["dataset"]]%in%x], simplify = FALSE)
              
              selectedFeatures <- mapply(.doSelection, 
                                         measurements = assayTrain,
                                         modellingParams = params[names(assayTrain)],
                                         MoreArgs = list(outcomesTrain = classesTrain, 
                                                         crossValParams = CrossValParams(permutations = 1, folds = 5), ###### Where to get this from?
                                                         verbose = 0)
              )
              
              do.call("rbind", selectedFeatures[2,])
              #S4Vectors::DataFrame(dataset = rep(names(selectedFeatures[2,]), unlist(lapply(selectedFeatures[2,], length))), feature = unlist(selectedFeatures[2,]))
          })
