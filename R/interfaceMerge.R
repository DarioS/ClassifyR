setGeneric("mergeTrainInterface", function(measurements, classes, ...)
{
    standardGeneric("mergeTrainInterface")
})


setMethod("mergeTrainInterface", "DFrame",
          function(measurements,
                   classes,
                   params,
                   ...)
          {
              assayTrain = split(data.frame(t(as.data.frame(measurements))), mcols(measurements)[["dataset"]])
              assayTrain = lapply(assayTrain, t)
              assayTrain = lapply(assayTrain, DataFrame)
              
              # if(length(classifierParams)>1){
              #     
              #     selectedFeatures <- mapply(ClassifyR:::.doSelection, 
              #                                measurements = measurements, 
              #                                selectParams = lapply(assayParams, function(x)x$selectParams),
              #                                trainParams = lapply(assayParams, function(x)x$trainParams), 
              #                                predictParams = lapply(assayParams, function(x)x$predictParams), 
              #                                MoreArgs = list(classes = classes, 
              #                                                training = seq_len(ncol(measurements$clinical)), 
              #                                                featureSets = NULL, 
              #                                                metaFeatures = NULL, 
              #                                                verbose = 0))
              # }else{
                  selectedFeatures <- mapply(ClassifyR:::.doSelection, 
                                             measurements = assayTrain,
                                             MoreArgs = list(classes = classes, 
                                                             training = seq_len(nrow(assayTrain[[1]])),
                                                             crossValParams = generateCrossValParams(nRepeats = 1, nFolds = 5, nCores = 1, selectionOptimisation= "Resubstitution"),
                                                             modellingParams = params,
                                                             verbose = 0)
                                             )
              # }
              
              
              
              # if ("clinical" %in% names(measurements)) {
                  
              features <- DataFrame(dataset = rep(names(selectedFeatures[2,]), unlist(lapply(selectedFeatures[2,], length))), feature = unlist(selectedFeatures[2,]))
              fullTrain <- measurements[, mcols(measurements) %in% features]
              

              selectParams = ClassifyR::SelectParams(
                  ClassifyR::differentMeansRanking,
                  tuneParams = list(nFeatures = ncol(fullTrain), 
                                    performanceType = "Balanced Error")
              )
              
              params2 <- params
              params2@selectParams <- selectParams
              
              
              runTestOutput = ClassifyR::runTest(
                  fullTrain,
                  classes,
                  training = seq_len(nrow(fullTrain)),
                  testing = seq_len(nrow(fullTrain)),
                  modellingParams = params2
              )
              
              
              fullModel = runTestOutput@models[[1]]
              
              fullModel
              
          })
