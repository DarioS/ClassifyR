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
      sampleFold <- rep(1:2, c(round(ncol(expression) * (100 - percent) / 100), round(ncol(expression) * percent / 100)))
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
    if(bootMode == "fold")
    {
      if(class(results[[1]][[1]][[2]]) == "list") multipleVarieties <- TRUE else multipleVarieties <- FALSE
    } else {
      if(class(results[[1]][[2]]) == "list") multipleVarieties <- TRUE else multipleVarieties <- FALSE
    }
    
    resultsLists <- lapply(c(4, 3, 2, 1, 5), function(position)
    {
      resultList <- lapply(results, function(sample)
      {
        if(position == 3 || multipleVarieties == FALSE)
        {
          if(!position %in% c(1, 2, 5))
          {
            if(bootMode == "fold")
            {
              if(position == 3 || class(sample[[1]][[4]]) != "data.frame")
                                  # First fold predictions are a not data.frame.
                unlist(lapply(sample, function(fold) fold[[position]]))
              else # Predictions is a data.frame with labels and scores.
                do.call(rbind, lapply(sample, function(fold) fold[[position]]))
            } else { #bootMode is "split".
              sample[[position]]
            }
          } else {
            if(bootMode == "fold")
              lapply(sample, function(fold) fold[[position]])
            else
              sample[[position]]
          }
        } else {
          if(position == 4) # Predictions with varieties.
          {
            if(bootMode == "fold")
            {
              prediction <- lapply(1:length(sample[[1]][[4]]), function(variety)
                            {
                              prediction <- lapply(sample, function(fold) fold[[position]][[variety]])
                              if(class(prediction[[1]]) == "data.frame")
                                do.call(rbind, prediction)
                              else
                                unlist(prediction)
                            })
            } else { # is 'split'.
            lapply(1:length(sample[[1]][[4]]), function(variety)
              sample[[position]][[variety]])
            }
          } else { # Ranked, selected features and test samples, with varieties.
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
      
      if(multipleVarieties == TRUE && position != 3) # Rearrange lists by variety.
        resultList <- lapply(1:length(resultList[[1]]), function(variety)
                        lapply(resultList, "[[", variety))
      else
        resultList
    })
  } else { # leave k out.
    if(class(results[[1]][[2]]) == "list") multipleVarieties <- TRUE else multipleVarieties <- FALSE
    resultsLists <- lapply(c(4, 3, 2, 1, 5), function(position)
    {
      resultList <- lapply(results, function(sample) sampleResult <- sample[[position]])
      
      if(position == 3 || position == 4 && multipleVarieties == FALSE && class(resultList[[1]]) != "data.frame")
        resultList <- unlist(resultList)
      else if (position == 4 && multipleVarieties == FALSE) # Predictions are in a data.frame.
        resultList <- do.call(rbind, resultList)
      else if(multipleVarieties == TRUE && position == 4 && class(resultList[[1]]) != "data.frame")
        lapply(1:length(resultList[[1]]), function(variety) unlist(lapply(resultList, "[[", variety)))
      else if(multipleVarieties == TRUE && position == 4)
        lapply(1:length(resultList[[1]]), function(variety) do.call(rbind, lapply(resultList, "[[", variety)))
      else if(multipleVarieties == TRUE && position %in% c(1, 2, 5) && is.function(selectionParams@featureSelection))
        lapply(1:length(resultList[[1]]), function(variety) lapply(resultList, "[[", variety))
      else
        resultList
    })
  }    
  
  if(multipleVarieties == FALSE)
    predictionsListed <- list(resultsLists[[1]]) else predictionsListed <- resultsLists[[1]]
  
  predictionTables <- lapply(1:length(predictionsListed), function(variety)
  {
    if(validation == "bootstrap")
    {
      lapply(1:length(predictionsListed[[1]]), function(resample)
      {
        switch(class(predictionsListed[[variety]][[resample]]),
               factor = data.frame(sample = resultsLists[[2]][[resample]], label = predictionsListed[[variety]][[resample]]),
               numeric = data.frame(sample = resultsLists[[2]][[resample]], score = predictionsListed[[variety]][[resample]]),
               data.frame = data.frame(sample = resultsLists[[2]][[resample]],
                                       label = predictionsListed[[variety]][[resample]][, sapply(predictionsListed[[variety]][[resample]], class) == "factor"],
                                       score = predictionsListed[[variety]][[resample]][, sapply(predictionsListed[[variety]][[resample]], class) == "numeric"]))
                
      })
    } else {
      switch(class(predictionsListed[[variety]]),
             factor = data.frame(sample = resultsLists[[2]], label = predictionsListed[[variety]]),
             numeric = data.frame(sample = resultsLists[[2]], score = predictionsListed[[variety]]),
             data.frame = data.frame(sample = resultsLists[[2]],
                                     label = predictionsListed[[variety]][, sapply(predictionsListed[[variety]], class) == "factor"],
                                     score = predictionsListed[[variety]][, sapply(predictionsListed[[variety]], class) == "numeric"]))
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

  if(multipleVarieties == FALSE)
  {
    if(length(unlist(resultsLists[[5]])) == 0) # All tune values are NULL.
      resultsLists[[5]] <- list(NULL)
    ClassifyResult(datasetName, classificationName, sampleNames(expression), featureNames(expression),
                   resultsLists[[4]], resultsLists[[3]], predictionTables[[1]],
                   pData(expression)[, "class"], validationInfo, resultsLists[[5]])
  } else {
    classifyResults <- lapply(1:length(predictionTables), function(variety)
    {
      if(length(unlist(resultsLists[[5]][[variety]])) == 0) # All tune values are NULL.
        resultsLists[[5]][[variety]] <- list(NULL)
      ClassifyResult(datasetName, classificationName, sampleNames(expression), featureNames(expression),
                     resultsLists[[4]][[variety]], resultsLists[[3]][[variety]],
                     predictionTables[[variety]], pData(expression)[, "class"], validationInfo, resultsLists[[5]][[variety]])
    })
    if(validation == "bootstrap")
      names(classifyResults) <- names(results[[1]][[1]][[4]])
    else
      names(classifyResults) <- names(results[[1]][[4]])
    classifyResults
  }
})
