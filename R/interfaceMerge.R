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
              measurements = split(as.data.frame(t(as.data.frame(measurements))), mcols(measurements)[["dataset"]])
              assayTrain = measurements #%>% lapply(as.matrix)
              
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
                                             measurements = measurement,
                                             MoreArgs = list(classes = classes, 
                                                             selectParams = params$selectParams, 
                                                             trainParams = params$trainParams, 
                                                             predictParams = params$predictParams, 
                                                             training = seq_len(nrow(measurements[[1]])),
                                                             featureSets = NULL, 
                                                             metaFeatures = NULL, 
                                                             verbose = 0)
                                             )
              # }
              
              
              
              # if ("clinical" %in% names(measurements)) {
              
              fullTrain <- mapply(function(x,y){x[,y]}, measurements, selectedFeatures[2,], SIMPLIFY = FALSE) %>%
                  do.call("cbind",.) %>%
                  DataFrame()
              
              
              selectParams = ClassifyR::SelectParams(
                  ClassifyR::differentMeansRanking,
                  "Combine",
                  tuneParams = list(nFeatures = ncol(fullTrain), 
                                    performanceType = "Balanced Error")
              )
              
              
              runTestOutput = ClassifyR::runTest(
                  fullTrain,
                  classes,
                  training = rownames(fullTrain),
                  testing = rownames(fullTrain),
                  params = list(selectParams, params$trainParams, params$predictParams),
                  datasetName = "merge",
                  classificationName = "merge"
              )
              
              
              fullModel = runTestOutput@models[[1]]
              
              fullModel
              
          })
