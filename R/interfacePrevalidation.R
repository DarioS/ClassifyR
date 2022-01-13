extractPrevalidation = function(assayPreval){ #}, startingCol) {
    # seq_len(length(assayPreval)) |>
    #     lapply(
    #         function(i)
    #             assayPreval[[i]] |>
    #             tibble::column_to_rownames("sample") |>
    #             dplyr::rename_all(function(col)
    #                 paste0(names(assayPreval)[i], col))
    #     ) |>
    #     # Taking all prevalidation class vectors except the last one (needed for multiple classes)
    #     lapply(function(x)
    #         x[, c(startingCol:(ncol(x) - 1)), drop = FALSE] %>% tibble::rownames_to_column("row_names")) |>
    #     purrr::reduce(merge, by = "row_names") |>
    #     janitor::clean_names() |>
    #     tibble::column_to_rownames("row_names")
    
    
    use <- which(names(assayPreval)!="clinical")
    
    assayPreval <- sapply(assayPreval[use], function(x){
        if(!"sample"%in%colnames(x))x <- data.frame(sample = rownames(x), x)
        x[order(x$sample),]}, simplify = FALSE)
    
    vec <- sapply(assayPreval, function(x){
        x[,!colnames(x) %in% c("sample", "permutation", "fold", "class")][,-1]
        }, simplify = TRUE)
    rownames(vec) <- assayPreval[[1]]$sample
    vec
}

# Use to pull out the names of features in a ClassifyR model - XYZ: Could Ditch if we really wanted to
featurePuller = function(classifyObject) {
    if ("selectResult" %in% slotNames(classifyObject)) {
        features = classifyObject@selectResult@chosenFeatures
    } else{
        features = classifyObject@originalFeatures
    }
}

setClass("prevalModel", slots = list(fullModel = "list"))


setGeneric("prevalTrainInterface", function(measurements, classes, ...)
{
    standardGeneric("prevalTrainInterface")
})


setMethod("prevalTrainInterface", "DFrame",
          function(measurements,
                   classes,
                   params,
                   ...)
          {
              
              assayTrain <- sapply(unique(mcols(measurements)[["dataset"]]), function(x) measurements[,mcols(measurements)[["dataset"]]%in%x], simplify = FALSE)
              
              if(! "clinical" %in% names(assayTrain)) stop("Must have a dataset called `clinical`")
              
              CVparams <- CrossValParams(permutations = 1, folds = 10, parallelParams = SerialParam(RNGseed = .Random.seed[1]), tuneMode = "Resubstitution") 
              
              usePreval <- names(assayTrain)[names(assayTrain)!="clinical"]
              assayTests <- bpmapply(
                  ClassifyR::runTests,
                  measurements = assayTrain[usePreval],
                  modellingParams = params[usePreval],
                  MoreArgs = list(
                      classes = classes,
                      crossValParams = CVparams,
                      verbose = 0
                  ), 
                  BPPARAM = SerialParam(RNGseed = .Random.seed[1])) |>
                  sapply(function(x)x@predictions, simplify = FALSE)
              
              prevalidationTrain <- extractPrevalidation(assayTests)
              
              # selectedFeaturesClinical <- runTest(assayTrain[["clinical"]],
              #                                    classes = classes,
              #                                    training = seq_len(nrow(assayTrain[["clinical"]])),
              #                                    testing = seq_len(nrow(assayTrain[["clinical"]])),
              #                                    modellingParams = params[["clinical"]],
              #                                    crossValParams = CVparams,
              #                                    .iteration = 1,
              #                                    verbose = 0
              # )$selected[, "feature"]
             
              #fullTrain = cbind(assayTrain[["clinical"]][,selectedFeaturesClinical], prevalidationTrain[rownames(assayTrain[["clinical"]]), , drop = FALSE])
              
              fullTrain = cbind(assayTrain[["clinical"]], prevalidationTrain[rownames(assayTrain[["clinical"]]), , drop = FALSE])
              
              
              

              
              finalModParam <- params[["clinical"]]
              #finalModParam@selectParams <- NULL
              
              runTestOutput = ClassifyR::runTest(
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
              
              prevalidationModels =  mapply(
                  ClassifyR::runTest,
                  measurements = assayTrain[usePreval],
                  modellingParams = params[usePreval],
                  MoreArgs = list(
                      classes = classes,
                      training = seq_len(nrow(fullTrain)),
                      testing =  seq_len(nrow(fullTrain)),
                      crossValParams = CVparams,
                      .iteration = 1,
                      verbose = 0
                  )
              )
              
              # Make a class
              fullModel$prevalidationModels <- prevalidationModels["models",]
              names(fullModel$prevalidationModels) <- colnames(prevalidationModels)
              fullModel$modellingParams <- params
              fullModel$prevalFeatures <- prevalidationModels["selected",] 
              
              fullModel <- new("prevalModel", fullModel = list(fullModel))
              fullModel
              
          })





setGeneric("prevalPredictInterface", function(fullModel, test, ...)
{
    standardGeneric("prevalPredictInterface")
})


setMethod("prevalPredictInterface", c("prevalModel", "DFrame"),
          function(fullModel,
                   test,
                   ...,
                   returnType = "both",
                   verbose = 0
          )
          {
              fullModel <- fullModel@fullModel[[1]]
              assayTest <- sapply(unique(mcols(test)[["dataset"]]), function(x) test[,mcols(test)[["dataset"]]%in%x], simplify = FALSE)
              
              prevalidationModels <- fullModel$prevalidationModels
              modelPredictionFunctions <- fullModel$modellingParams
              
              prevalidationPredict <- sapply(names(prevalidationModels), function(x){
                  predictParams <- modelPredictionFunctions[[x]]@predictParams
                  paramList <- list(prevalidationModels[[x]], assayTest[[x]])
                  if(length(predictParams@otherParams) > 0) paramList <- c(paramList, predictParams@otherParams)
                  paramList <- c(paramList, verbose = 0)
                  prediction <- do.call(predictParams@predictor, paramList)
                  prediction}, simplify = FALSE) |>
                  extractPrevalidation()
              
 
              fullTest = cbind(assayTest[["clinical"]], prevalidationPredict[rownames(assayTest[["clinical"]]), , drop = FALSE])
              
              
              predictParams <- modelPredictionFunctions[["clinical"]]@predictParams
              paramList <- list(fullModel,  fullTest)
              if(length(predictParams@otherParams) > 0) paramList <- c(paramList, predictParams@otherParams)
              paramList <- c(paramList, verbose = 0)
              finalPredictions <- do.call(predictParams@predictor, paramList)

              finalPredictions
          })
