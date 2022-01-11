setGeneric("selectMulti", function(measurements, classes, params, ...)
    standardGeneric("selectMulti"))

# DataFrame of numeric measurements, likely created by runTests or runTest.
setMethod("selectMulti", "DataFrame",
          function(measurements, classes, params, verbose = 0)
          {
              
              assayTrain <- split(data.frame(t(as.data.frame(measurements))), mcols(measurements)[["dataset"]])
              assayTrain <- lapply(assayTrain, t)
              assayTrain <- lapply(assayTrain, DataFrame)
              
              selectedFeatures <- mapply(ClassifyR:::.doSelection, 
                                         measurements = assayTrain,
                                         modellingParams = params[names(assayTrain)],
                                         MoreArgs = list(classes = classes, 
                                                         training = seq_len(nrow(assayTrain[[1]])),
                                                         crossValParams = crossValParams,
                                                         verbose = 0)
              )
              
              DataFrame(dataset = rep(names(selectedFeatures[2,]), unlist(lapply(selectedFeatures[2,], length))), feature = unlist(selectedFeatures[2,]))
          })
