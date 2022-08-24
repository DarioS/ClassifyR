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
    #         x[, c(startingCol:(ncol(x) - 1)), drop = FALSE] |> tibble::rownames_to_column("row_names")) |>
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

prevalTrainInterface <- function(measurements, classes, params, ...)
          {
              
              ###
              # Splitting measurements into a list of each of the assays
              ###
              assayTrain <- sapply(unique(S4Vectors::mcols(measurements)[["assay"]]), function(assay) measurements[, S4Vectors::mcols(measurements)[["assay"]] %in% assay], simplify = FALSE)
              
              if(!"clinical" %in% names(assayTrain)) stop("Must have an assay called \"clinical\"")
              
              # Create generic crossValParams to use for my prevalidation... should add this as a input parameter
              CVparams <- CrossValParams(permutations = 1, folds = 10, parallelParams = SerialParam(RNGseed = .Random.seed[1]), tuneMode = "Resubstitution") 
              
              ###
              # Fit a classification model for each non-clinical datasets, pulling models from "params"
              ###
              usePreval <- names(assayTrain)[names(assayTrain) != "clinical"]
              assayTests <- bpmapply(
                  runTests,
                  measurements = assayTrain[usePreval],
                  modellingParams = params[usePreval],
                  MoreArgs = list(
                      outcome = classes,
                      crossValParams = CVparams,
                      verbose = 0
                  ), 
                  BPPARAM = SerialParam(RNGseed = .Random.seed[1])) |>
                  sapply(function(result) result@predictions, simplify = FALSE)
              
              ###
              # Pull-out prevalidated vectors ie. the predictions on each of the test folds.
              ###
              prevalidationTrain <- extractPrevalidation(assayTests)
              
              # Feature select on clinical data before binding
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
              
              prevalidationTrain <- S4Vectors::DataFrame(prevalidationTrain)
              mcols(prevalidationTrain)$assay = "prevalidation"
              mcols(prevalidationTrain)$feature = colnames(prevalidationTrain)
              
              
              ###
              # Bind the prevalidated data to the clinical data
              ###
              fullTrain = cbind(assayTrain[["clinical"]], prevalidationTrain[rownames(assayTrain[["clinical"]]), , drop = FALSE])
              
              
              

              # Pull out clinical data
              finalModParam <- params[["clinical"]]
              #finalModParam@selectParams <- NULL
              
              # Fit classification model (from clinical in params)
              runTestOutput = runTest(
                  measurementsTrain = fullTrain,
                  outcomeTrain = classes,
                  measurementsTest = fullTrain,
                  outcomeTest = classes,
                  modellingParams = finalModParam,
                  crossValParams = CVparams,
                  .iteration = 1,
                  verbose = 0
                  )
              
              
              # Extract the classification model from runTest output
              fullModel = runTestOutput$models
              fullModel$fullFeatures = colnames(fullTrain)
              
              # Fit models with each datatype for use in prevalidated prediction later..
              prevalidationModels =  mapply(
                  runTest,
                  measurementsTrain = assayTrain,
                  measurementsTest = assayTrain,               
                  modellingParams = params,
                  MoreArgs = list(
                      crossValParams = CVparams,
                      outcomeTrain = classes,
                      outcomeTest = classes,
                      .iteration = 1,
                      verbose = 0
                  )
              )
              
              # Add prevalidated models and classification params for each datatype to the fullModel object
              fullModel$prevalidationModels <- prevalidationModels["models",]
              names(fullModel$prevalidationModels) <- colnames(prevalidationModels)
              fullModel$modellingParams <- params
              fullModel$prevalFeatures <- prevalidationModels["selected",] 
              
              fullModel <- new("prevalModel", fullModel = list(fullModel))
              fullModel
}

prevalPredictInterface <- function(fullModel, test, ..., returnType = "both", verbose = 0)
          {
              fullModel <- fullModel@fullModel[[1]]
              assayTest <- sapply(unique(mcols(test)[["assay"]]), function(assay) test[, mcols(test)[["assay"]] %in% assay], simplify = FALSE)
              
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
              
              prevalidationPredict <- S4Vectors::DataFrame(prevalidationPredict)
              mcols(prevalidationPredict)$assay = "prevalidation"
              mcols(prevalidationPredict)$feature = colnames(prevalidationPredict)
              
              fullTest = cbind(assayTest[["clinical"]], prevalidationPredict[rownames(assayTest[["clinical"]]), , drop = FALSE])
              
              
              predictParams <- modelPredictionFunctions[["clinical"]]@predictParams
              paramList <- list(fullModel,  fullTest)
              if(length(predictParams@otherParams) > 0) paramList <- c(paramList, predictParams@otherParams)
              paramList <- c(paramList, verbose = 0)
              finalPredictions <- do.call(predictParams@predictor, paramList)

              finalPredictions
          }