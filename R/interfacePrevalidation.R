extractPrevalidation = function(assayPreval, startingCol) {
    seq_len(length(assayPreval)) %>%
        lapply(
            function(i)
                assayPreval[[i]] %>%
                tibble::column_to_rownames("sample") %>%
                dplyr::rename_all(function(col)
                    paste0(names(assayPreval)[i], col))
        ) %>%
        # Taking all prevalidation class vectors except the last one (needed for multiple classes)
        lapply(function(x)
            x[, c(startingCol:(ncol(x) - 1)), drop = FALSE] %>% tibble::rownames_to_column("row_names")) %>%
        purrr::reduce(merge, by = "row_names") %>%
        janitor::clean_names() %>%
        tibble::column_to_rownames("row_names")
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
                   foldParams,
                   performanceType,
                   seed,
                   ...)
          {
       
              assayTrain <- sapply(unique(mcols(measurements)[["dataset"]]), function(x) measurements[,mcols(measurements)[["dataset"]]%in%x], simplify = FALSE)
              
              
              assayParams = mapply(
                  FUN = paramGenerator,
                  assayTrain,
                  classifier = classifierParams[names(classifierParams) != "clinical"],
                  featureSelectionMethod =  featureSelectionParams[names(featureSelectionParams) != "clinical"],
                  nFeatures = nFeatureParams[names(nFeatureParams) != "clinical"],
                  MoreArgs = list(performanceType = performanceType,
                                  params = NULL),
                  SIMPLIFY = FALSE
              )
              
              
              assayPreval = mapply(
                  ClassifyR::runTests,
                  assayTrain,
                  classificationName = paste("prevalidation", names(assayTrain)),
                  datasetName = names(assayTrain),
                  params = assayParams,
                  folds = foldParams[names(foldParams) != "clinical"],
                  MoreArgs = list(
                      classes = classes,
                      permutations = 1,
                      seed = seed
                  )
              ) %>% lapply(function(x)
                  x@predictions) %>% purrr::flatten()
              
              prevalidationTrain = extractPrevalidation(assayPreval, startingCol = 2)
              
              paramsClinical = paramGenerator(
                  measurements$clinical,
                  classifier = classifierParams$clinical,
                  featureSelectionMethod =  featureSelectionParams$clinical,
                  nFeatures = nFeatureParams$clinical,
                  performanceType = performanceType,
                  params = NULL
              )
              
              selectedFeaturesClinical = runTest(as.matrix(measurements$clinical),
                                                 classes = classes,
                                                 training = colnames(measurements$clinical),
                                                 testing =  colnames(measurements$clinical),
                                                 params = paramsClinical,
                                                 classificationName = "clinical",
                                                 datasetName = "clinical",
              )@selectResult@rankedFeatures[[1]]
              
              # if ("clinical" %in% names(measurements)) {
              
              fullTrain = data.frame(t(measurements[["clinical"]]))[,selectedFeaturesClinical] %>%
                  cbind(prevalidationTrain[colnames(measurements[["clinical"]]), , FALSE])
              
              # } else{
              #     fullTrain = cbind(prevalidationTrain)
              # }
              
              
              # XYZ - should the feature selection be here
              # paramsFull = paramGenerator(
              #     t(fullTrain),
              #     classifier = classifierParams$clinical,
              #     featureSelectionMethod = featureSelectionParams$clinical,
              #     nFeatures = nFeatureParams$clinical,
              #     performanceType = performanceType
              # )
              
              paramsFull = paramGenerator(
                  DataFrame(fullTrain),
                  classifier = classifierParams$clinical,
                  featureSelectionMethod = NULL,
                  nFeatures = NULL,
                  performanceType = performanceType
              )
              
              
              runTestOutput = ClassifyR::runTest(
                  DataFrame(fullTrain),
                  classes,
                  training = rownames(fullTrain),
                  testing = rownames(fullTrain),
                  params = paramsFull,
                  datasetName = "prevalidation",
                  classificationName = "prevalidation"
              )
              
              
              fullModel = runTestOutput@models[[1]]
              # XYZ - feature selection happening even when it theoretically shouldn't
              # fullModel$fullFeatures = runTestOutput %>% featurePuller() %>% unlist()
              fullModel$fullFeatures = runTestOutput@originalFeatures
              
              prevalidationModels =  mapply(
                  ClassifyR::runTest,
                  measurements = assayTrain,
                  classificationName = paste("prevalidation", names(assayTrain)),
                  datasetName = names(assayTrain),
                  params = assayParams,
                  MoreArgs = list(
                      classes = classes,
                      training = colnames(assayTrain[[1]]),
                      testing =  colnames(assayTrain[[1]])
                  )
              )
              
              # Make a class
              fullModel$prevalidationModels = prevalidationModels %>% lapply(function(x)
                  x@models) %>% purrr::flatten()
              fullModel$fullParams = paramsFull$predictParams
              fullModel$prevalParams = assayParams
              fullModel$prevalFeatures = prevalidationModels %>% lapply(featurePuller) %>% purrr::flatten()
              fullModel$assayParams = assayParams
              
              fullModel = new("prevalModel", fullModel = list(fullModel))
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
              fullModel = fullModel@fullModel[[1]]
              test = split(as.data.frame(t(as.data.frame(test))), mcols(test)[["dataset"]])
              prevalidationModels = fullModel$prevalidationModels
              prevalListModels = prevalidationModels %>% lapply(function(x)
                  list(x))
              modelPredictionFunctions = fullModel$assayParams %>% lapply(function(x)
                  x$predictParams@predictor)
              
              
              testAssays = test[names(test) != "clinical"] %>% lapply(function(x)
                  x %>% t() %>%  data.frame())
              
              testAssays = testAssays %>% lapply(function(x)
                  x %>% DataFrame)
              
              predictedPrevalidation = mapply(
                  lapply,
                  X = prevalListModels ,
                  modelPredictionFunctions,
                  test = testAssays,
                  MoreArgs = list(returnType = "both", verbose = verbose),
                  SIMPLIFY = FALSE
              ) %>% purrr::flatten() %>%
                  lapply(function(x)
                      data.frame(x) %>% tibble::rownames_to_column("sample")) %>%
                  extractPrevalidation(2)
              
              # if ("clinical" %in% names(test)) {
              fullTest = data.frame(t(test[["clinical"]])) %>%
                  cbind(predictedPrevalidation[colnames(test[["clinical"]]), , FALSE])
              # } else{
              #     fullTest = predictedPrevalidation
              # }
              
              
              finalPredictions = fullModel$fullParams@predictor(fullModel,
                                                                DataFrame(fullTest),
                                                                returnType = returnType,
                                                                verbose = verbose)
              finalPredictions
          })
