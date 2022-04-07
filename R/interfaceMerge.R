setGeneric("mergeTrainInterface", function(measurements, classes, params, ...)
{
    standardGeneric("mergeTrainInterface")
})


setMethod("mergeTrainInterface", "DFrame",
          function(measurements,
                   classes,
                   params,
                   ...)
          {
              fullTrain <- measurements
              
              params2 <- params
              params2@selectParams <- NULL
              
              
              runTestOutput = runTest(
                  fullTrain,
                  classes,
                  training = seq_len(nrow(fullTrain)),
                  testing = seq_len(nrow(fullTrain)),
                  modellingParams = params2
              )
              
              
              fullModel = runTestOutput@models[[1]]
              
              fullModel
              
          })
