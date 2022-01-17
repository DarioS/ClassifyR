setGeneric("selectMulti", function(measurements, classes, params, ...)
    standardGeneric("selectMulti"))

# DataFrame of numeric measurements, likely created by runTests or runTest.
setMethod("selectMulti", "DataFrame",
          function(measurements, classes, params, verbose = 0)
          {
              
              assayTrain <- sapply(unique(mcols(measurements)[["dataset"]]), function(x) measurements[,mcols(measurements)[["dataset"]]%in%x], simplify = FALSE)
              
              selectedFeatures <- mapply(.doSelection, 
                                         measurements = assayTrain,
                                         modellingParams = params[names(assayTrain)],
                                         MoreArgs = list(classes = classes, 
                                                         training = seq_len(nrow(assayTrain[[1]])),
                                                         crossValParams = CrossValParams(permutations = 1, folds = 5), ###### Where to get this from?
                                                         verbose = 0)
              )
              
              do.call("rbind", selectedFeatures[2,])
              #DataFrame(dataset = rep(names(selectedFeatures[2,]), unlist(lapply(selectedFeatures[2,], length))), feature = unlist(selectedFeatures[2,]))
          })
