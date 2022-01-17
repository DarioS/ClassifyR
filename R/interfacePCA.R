
setClass("pcaModel", slots = list(fullModel = "list"))


setGeneric("pcaTrainInterface", function(measurements, classes, ...)
{
    standardGeneric("pcaTrainInterface")
})


setMethod("pcaTrainInterface", "DFrame",
          function(measurements,
                   classes,
                   params,
                   nFeatures,
                   ...)
          {
              
              assayTrain <- sapply(unique(mcols(measurements)[["dataset"]]), function(x) measurements[,mcols(measurements)[["dataset"]]%in%x], simplify = FALSE)
              
              if(! "clinical" %in% names(assayTrain)) stop("Must have a dataset called `clinical`")
              
              CVparams <- CrossValParams(permutations = 1, folds = 10, parallelParams = SerialParam(RNGseed = .Random.seed[1]), tuneMode = "Resubstitution") 
              
              usePCA<- names(assayTrain)[names(assayTrain)!="clinical"]
              assayPCA <- sapply(assayTrain[usePCA], function(assay){
                  pca <- prcomp(assay)
              }, simplify = FALSE)
              
              
              pcaVar <- mapply(function(pca, n){
                  data.frame(pca$x[,seq_len(n)])
              }, assayPCA, nFeatures[names(assayPCA)], SIMPLIFY = FALSE)
              
              pcaVar <- do.call("cbind", pcaVar)
              
              pcaVar <- DataFrame(pcaVar)
              mcols(pcaVar)$dataset = "pca"
              mcols(pcaVar)$feature = colnames(pcaVar)
            
              fullTrain = cbind(assayTrain[["clinical"]], pcaVar)
              
              
              
              finalModParam <- params[["clinical"]]
              
              runTestOutput = runTest(
                  fullTrain,
                  classes = classes,
                  training = seq_len(nrow(fullTrain)),
                  testing = seq_len(nrow(fullTrain)),
                  modellingParams = finalModParam,
                  crossValParams = CVparams,
                  .iteration = 1,
                  verbose = 0
                  )
              
              
              fullModel = runTestOutput$models
              fullModel$fullFeatures = colnames(fullTrain)
              
            
              # Make a class
              fullModel$pcaModels <- assayPCA
              fullModel$nFeatures <- nFeatures
              fullModel$modellingParams <- finalModParam
              fullModel <- new("pcaModel", fullModel = list(fullModel))
              fullModel
              
          })





setGeneric("pcaPredictInterface", function(fullModel, test, ...)
{
    standardGeneric("pcaPredictInterface")
})


setMethod("pcaPredictInterface", c("pcaModel", "DFrame"),
          function(fullModel,
                   test,
                   ...,
                   returnType = "both",
                   verbose = 0
          )
          {
              fullModel <- fullModel@fullModel[[1]]
              assayTest <- sapply(unique(mcols(test)[["dataset"]]), function(x) test[,mcols(test)[["dataset"]]%in%x], simplify = FALSE)
              
              pcaModels <- fullModel$pcaModels
              nFeatures <- fullModel$nFeatures
              
              pcaVar <- mapply(function(pca, assay, n){
                 data.frame(predict(pca, assay))[,seq_len(n)]
              }, pcaModels, assayTest[names(pcaModels)], nFeatures[names(pcaModels)], SIMPLIFY = FALSE)
              
              pcaVar <- do.call(cbind, pcaVar)
              
              pcaVar <- DataFrame(pcaVar)
              mcols(pcaVar)$dataset = "pca"
              mcols(pcaVar)$feature = colnames(pcaVar)
              
              fullTest = cbind(assayTest[["clinical"]], pcaVar)

              
              modelPredictionFunctions <- fullModel$modellingParams
              predictParams <- modelPredictionFunctions@predictParams
              paramList <- list(fullModel,  fullTest)
              if(length(predictParams@otherParams) > 0) paramList <- c(paramList, predictParams@otherParams)
              paramList <- c(paramList, verbose = 0)
              finalPredictions <- do.call(predictParams@predictor, paramList)

              finalPredictions
          })
