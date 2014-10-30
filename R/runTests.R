setGeneric("runTests", function(expression, ...)
           {standardGeneric("runTests")})

setMethod("runTests", c("matrix"),
  function(expression, classes, ...)
{
    groupsTable <- data.frame(class = classes)
    rownames(groupsTable) <- colnames(expression)
    exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
    runTests(exprSet, ...)
})

setMethod("runTests", c("ExpressionSet"),
    function(expression, datasetName, classificationName, validation = c("bootstrap", "leaveOut"),
             bootMode = c("fold", "split"),
             resamples = 100, percent = 25, folds = 5, leave = 2,
             seed, parallelParams = bpparam(),
             params = list(SelectionParams(), TrainParams(), PredictParams()),
             verbose = 1)
{
  validation <- match.arg(validation)
  bootMode <- match.arg(bootMode)
  if(!missing(seed)) set.seed(seed)

  stagesParamClasses <- sapply(params, class)
  if(match("TrainParams", stagesParamClasses) > match("PredictParams", stagesParamClasses))
    stop("\"testing\" must not be before \"training\" in 'params'.")
  
  if(validation == "bootstrap")
  {
    bootstrap <- lapply(1:resamples, function(bootstrap) sample(ncol(expression), replace = TRUE))
    if(bootMode == "fold")
    {
      sampleFold <- rep(1:folds, length.out = ncol(expression))
      samplesFolds <- lapply(bootstrap, function(sample) split(sample, sampleFold))
    } else {
      sampleFold <- rep(1:2, c(ncol(expression) * (100 - percent) / 100, ncol(expression) * percent / 100))
      samplesFolds <- lapply(bootstrap, function(sample) split(sample, sampleFold))
    }

    results <- bpmapply(function(sampleFolds, sampleNumber)
    {
      if(verbose >= 1 && sampleNumber %% 10 == 0)
        message("Processing sample set ", sampleNumber, '.')
      if(bootMode == "fold")
      {
        lapply(1:length(sampleFolds), function(foldIndex)
        {
          runTest(expression, training = unlist(sampleFolds[-foldIndex]),
                  testing = sampleFolds[[foldIndex]], params = params, verbose = verbose)
        })
      } else { # Split mode.
        runTest(expression, training = sampleFolds[[1]],
                testing = sampleFolds[[2]], params = params, verbose = verbose)
      }
    }, samplesFolds, as.list(1:length(samplesFolds)), BPPARAM = parallelParams, SIMPLIFY = FALSE)
  } else # leave k out.
  {
    testSamples <- as.data.frame(combn(ncol(expression), leave))
    trainingSamples <- lapply(testSamples, function(sample) setdiff(1:ncol(expression), sample))
    results <- bpmapply(function(trainingSample, testSample, sampleNumber)
    {
      if(verbose >= 1 && sampleNumber %% 10 == 0)
        message("Processing sample set ", sampleNumber, '.')
      runTest(expression, training = trainingSample, testing = testSample,
              params = params, verbose = verbose)
    }, trainingSamples, testSamples, (1:length(trainingSamples)),
    BPPARAM = parallelParams, SIMPLIFY = FALSE)
  }

  selectionParams <- params[[match("SelectionParams", stagesParamClasses)]]
  predictParams <- params[[match("PredictParams", stagesParamClasses)]]
  if(validation == "bootstrap")
  {
    resultsLists <- lapply(c(4, 3, 2, 1), function(position)
    {
      resultList <- lapply(results, function(sample)
      {
        
        if(position == 3 || predictParams@multipleResults == FALSE)
        {
          if(!position %in% c(1, 2))
          {
            if(bootMode == "fold")
              unlist(lapply(sample, function(fold) fold[[position]]))
            else
              unlist(sample[[position]])
          } else {
            if(bootMode == "fold")
              lapply(sample, function(fold) fold[[position]])
            else
              sample[[position]]
          }
        } else {
          if(position == 4)
          {
            if(bootMode == "fold")
            lapply(1:length(sample[[1]][[4]]), function(variety)
              unlist(lapply(sample, function(fold) fold[[position]][[variety]])))
            else
            lapply(1:length(sample[[1]][[4]]), function(variety)
              unlist(factor(sample[[position]][[variety]])))
          } else {
            if(is.function(selectionParams@featureSelection))
            {
              if(bootMode == "fold")
                lapply(1:length(sample[[1]][[4]]), function(variety)
                  lapply(sample, function(fold) fold[[position]][[variety]]))
              else
                lapply(1:length(sample[[1]][[4]]), function(variety)
                  sample[[position]][[variety]])
            } else {
              if(bootMode == "fold")
                lapply(sample, function(fold) fold[[position]])
              else
                sample[["position"]]
            }
          }
        }
      })
      
      if(predictParams@multipleResults == TRUE && position %in% c(1, 2, 4))
        resultList <- lapply(1:length(resultList[[1]]), function(variety)
                        lapply(resultList, "[[", variety))
      else
        resultList
    })
  } else { # leave k out.
    resultsLists <- lapply(c(4, 3, 2, 1), function(position)
    {
      resultList <- lapply(results, function(sample) sampleResult <- sample[[position]])
      
      if(position == 3 || (predictParams@multipleResults == FALSE && position == 4))
        resultList <- unlist(resultList)
      else if(predictParams@multipleResults == TRUE && position == 4)
        lapply(1:length(resultList[[1]]), function(variety) unlist(lapply(resultList, "[[", variety)))
      else if(predictParams@multipleResults == TRUE && position %in% c(1, 2) && is.function(selectionParams@featureSelection))
        lapply(1:length(resultList[[1]]), function(variety) lapply(resultList, "[[", variety))
      else
        resultList
    })
  }    
  
  if(predictParams@multipleResults == FALSE)
    predictionsListed <- list(resultsLists[[1]]) else predictionsListed <- resultsLists[[1]]
  
  predictionTables <- lapply(1:length(predictionsListed), function(variety)
  {
    if(validation == "bootstrap")
    {
      lapply(1:length(predictionsListed[[1]]), function(resample)
             data.frame(sample = resultsLists[[2]][[resample]], predicted = predictionsListed[[variety]][[resample]]))        
    } else {
      list(data.frame(sample = resultsLists[[2]], predicted = predictionsListed[[variety]]))
    }
  })
  
  if(validation == "bootstrap")
  {
    if(bootMode == "fold")
      validationInfo <- list("fold", resamples, folds)
    else
      validationInfo <- list("split", resamples, percent)
  } else {
    validationInfo <- list("leave", leave)
  }

  if(predictParams@multipleResults == FALSE)
  {
    ClassifyResult(datasetName, classificationName, sampleNames(expression), featureNames(expression),
                   resultsLists[[4]], resultsLists[[3]], predictionTables[[1]],
                   pData(expression)[, "class"], validationInfo)
  } else {
    classifyResults <- lapply(1:length(predictionTables), function(variety)
    {
      ClassifyResult(datasetName, classificationName, sampleNames(expression), featureNames(expression),
                     resultsLists[[4]][[variety]], resultsLists[[3]][[variety]],
                     predictionTables[[variety]], pData(expression)[, "class"], validationInfo)
    })
    if(validation == "bootstrap")
      names(classifyResults) <- names(results[[1]][[1]][[4]])
    else
      names(classifyResults) <- names(results[[1]][[4]])
    classifyResults
  }
})
