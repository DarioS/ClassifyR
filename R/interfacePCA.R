
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
              
              ###
              # Splitting measurements into a list of each of the datasets
              ###
              assayTrain <- sapply(unique(mcols(measurements)[["dataset"]]), function(x) measurements[,mcols(measurements)[["dataset"]]%in%x], simplify = FALSE)
              
              if(! "clinical" %in% names(assayTrain)) stop("Must have a dataset called `clinical`")
              
              # Create generic crossValParams just to get things working, might be used for optimising features in runTest later???
              CVparams <- CrossValParams(permutations = 1, folds = 10, parallelParams = SerialParam(RNGseed = .Random.seed[1]), tuneMode = "Resubstitution") 
              
              ###
              # Run PCA for all datasets except clinical
              ###
              usePCA<- names(assayTrain)[names(assayTrain)!="clinical"]
              assayPCA <- sapply(assayTrain[usePCA], function(assay){
                  pca <- prcomp(assay)
              }, simplify = FALSE)
              
              
              ###
              # Pull out the top n (nFeatures) PCA loadings
              ###
              
              pcaVar <- mapply(function(pca, n){
                  data.frame(pca$x[,seq_len(n)])
              }, assayPCA, nFeatures[names(assayPCA)], SIMPLIFY = FALSE)
              
              pcaVar <- do.call("cbind", pcaVar)
              
              ###
              # Bind the PCA components to the clinical data
              ###
              
              pcaVar <- DataFrame(pcaVar)
              mcols(pcaVar)$dataset = "pca"
              mcols(pcaVar)$feature = colnames(pcaVar)
            
              fullTrain = cbind(assayTrain[["clinical"]], pcaVar)
              
              
              ###
              # Pull out the modellingParams for the clinical data.
              ###
              
              finalModParam <- params[["clinical"]]
              
              
              ###
              # Fit model with runTest
              ### 
              
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
              
              ###
              # Pull out fitted model from runTest
              ###
              
              fullModel = runTestOutput$models
              fullModel$fullFeatures = colnames(fullTrain)
              
            
              # Make a class
              fullModel$pcaModels <- assayPCA # Add PCA "models" to final model.
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
              # Pull out my classification model
              fullModel <- fullModel@fullModel[[1]]
              
              #Split my test data into a list of the different datasets
              assayTest <- sapply(unique(mcols(test)[["dataset"]]), function(x) test[,mcols(test)[["dataset"]]%in%x], simplify = FALSE)
              
              # Pull out my PCA models
              pcaModels <- fullModel$pcaModels
              nFeatures <- fullModel$nFeatures
              
              # Project test data into trained PCA space
              pcaVar <- mapply(function(pca, assay, n){
                 data.frame(predict(pca, assay))[,seq_len(n)]
              }, pcaModels, assayTest[names(pcaModels)], nFeatures[names(pcaModels)], SIMPLIFY = FALSE)
              
              pcaVar <- do.call(cbind, pcaVar)
              
              pcaVar <- DataFrame(pcaVar)
              mcols(pcaVar)$dataset = "pca"
              mcols(pcaVar)$feature = colnames(pcaVar)
              
              # Merge my PCA stuff with my clinical data
              fullTest = cbind(assayTest[["clinical"]], pcaVar)

              # Pull out my predictParams for clinical
              modelPredictionFunctions <- fullModel$modellingParams
              predictParams <- modelPredictionFunctions@predictParams
              paramList <- list(fullModel,  fullTest)
              if(length(predictParams@otherParams) > 0) paramList <- c(paramList, predictParams@otherParams)
              paramList <- c(paramList, verbose = 0)
              
              # Predict with my predictor
              finalPredictions <- do.call(predictParams@predictor, paramList)

              finalPredictions
          })
