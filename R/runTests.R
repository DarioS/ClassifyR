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
    function(expression, datasetName, classificationName, validation = c("bootstrap", "leaveOut", "fold"),
             bootMode = c("fold", "split"),
             resamples = 100, percent = 25, folds = 5, leave = 2,
             seed, parallelParams = bpparam(),
             params = list(SelectParams(), TrainParams(), PredictParams()),
             verbose = 1)
{
  validation <- match.arg(validation)
  bootMode <- match.arg(bootMode)
  if(!missing(seed)) set.seed(seed)
  resultTypes <- c("ranked", "selected", "testSet", "predictions", "tune")

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
        lapply(1:folds, function(foldIndex)
        {
          runTest(expression, training = unlist(sampleFolds[-foldIndex]),
                  testing = sampleFolds[[foldIndex]], params = params, verbose = verbose,
                  .iteration = c(sampleNumber, foldIndex))
        })
      } else { # Split mode.
        runTest(expression, training = sampleFolds[[1]],
                testing = sampleFolds[[2]], params = params, verbose = verbose, .iteration = sampleNumber)
      }
    }, samplesFolds, as.list(1:resamples), BPPARAM = parallelParams, SIMPLIFY = FALSE)
  } else if(validation == "leave") # leave k out.
  {
    testSamples <- as.data.frame(combn(ncol(expression), leave))
    trainingSamples <- lapply(testSamples, function(sample) setdiff(1:ncol(expression), sample))
    results <- bpmapply(function(trainingSample, testSample, sampleNumber)
    {
      if(verbose >= 1 && sampleNumber %% 10 == 0)
        message("Processing sample set ", sampleNumber, '.')
      runTest(expression, training = trainingSample, testing = testSample,
              params = params, verbose = verbose, .iteration = sampleNumber)
    }, trainingSamples, testSamples, (1:length(trainingSamples)),
    BPPARAM = parallelParams, SIMPLIFY = FALSE)
  } else { # Unresampled, ordinary k-fold cross-validation.
    sampleFold <- rep(1:folds, length.out = ncol(expression))
    samplesFolds <- split(1:ncol(expression), sampleFold)
    
    if(verbose >= 1)
      message("Processing ", folds, "-fold cross-validation.")
    results <- bplapply(1:folds, function(foldIndex)
    {
      runTest(expression, training = unlist(samplesFolds[-foldIndex]),
              testing = samplesFolds[[foldIndex]], params = params, verbose = verbose,
              .iteration = foldIndex)
    }, BPPARAM = parallelParams)
  }
  
  if(validation == "bootstrap" && bootMode == "split" || validation %in% c("leaveOut", "fold"))
  {
    resultErrors <- sapply(results, function(result) is.character(result))
    if(sum(resultErrors) == length(results))
    {
      message("Error: All cross-validations had an error.")
      return(results)
    } else if(sum(resultErrors) != 0) # Filter out cross-validations resulting in error.
    {
      results <- results[!resultErrors]
    }
  } else { # Result has nested lists, because of bootstrap resampling with folding.
    resultErrors <- lapply(results, function(resample) lapply(resample, is.character))
    if(sum(unlist(resultErrors)) == resamples * folds)
    {
      message("Error: All cross-validations had an error.")
      return(results)
    } else if(sum(unlist(resultErrors)) != 0) # Filter out error cross-validations.
    {
      results <- results[sapply(results, function(resample) !all(sapply(resample, is.character)))]
      results <- lapply(results, function(resample) resample[!sapply(resample, is.character)])
    }
  }

  selectParams <- params[[match("SelectParams", stagesParamClasses)]]
  predictParams <- params[[match("PredictParams", stagesParamClasses)]]
  if(validation == "bootstrap")
  {
    if(bootMode == "fold")
    {
      if(class(results[[1]][[1]][["predictions"]]) == "list") 
      {
        multipleVarieties <- TRUE
        varietyNames <- names(results[[1]][[1]][["predictions"]])
      } else {
        multipleVarieties <- FALSE
        varietyNames <- "None"
      }
    } else {
      if(class(results[[1]][["predictions"]]) == "list")
      {
        multipleVarieties <- TRUE
        varietyNames <- names(results[[1]][["predictions"]])
      } else {
        multipleVarieties <- FALSE
        varietyNames <- "None"
      }
    }

    resultsByVariety <- lapply(varietyNames, function(varietyName)
    {
      resultsLists <- lapply(resultTypes, function(resultType)
      {
        resultList <- lapply(results, function(sample)
        {
          if(resultType == "testSet" || multipleVarieties == FALSE)
          {
            if(resultType %in% c("testSet", "predictions"))
            {
              if(bootMode == "fold")
              {
                if(resultType == "testSet" || class(sample[[1]][["predictions"]]) != "data.frame")
                                              # First fold predictions are a not data.frame.
                  unlist(lapply(sample, function(fold) fold[[resultType]]))
                else # Predictions is a data.frame with labels and scores.
                  do.call(rbind, lapply(sample, function(fold) fold[[resultType]]))
              } else { # bootMode is "split".
                sample[[resultType]]
              }
            } else {
              if(bootMode == "fold")
                lapply(sample, function(fold) fold[[resultType]])
              else
                sample[[resultType]]
            }
          } else { # There are multiple varieties.
            if(resultType == "predictions") # Predictions with varieties.
            {
              if(bootMode == "fold")
              {
                prediction <- lapply(sample, function(fold) fold[[resultType]][[varietyName]])
                if(class(prediction[[1]]) == "data.frame")
                  do.call(rbind, prediction)
                else
                  unlist(prediction)
              } else { # bootMode is "split".
                sample[[resultType]][[varietyName]]
              }
            } else { # Ranked, selected features and test samples, with varieties.
              if(bootMode == "fold")
              {
                if(is.null(sample[[1]][[resultType]]))
                  NULL
                else
                  lapply(sample, function(fold) fold[[resultType]][[varietyName]])
              } else {
                if(is.null(sample[[resultType]]))
                  NULL
                else
                  sample[[resultType]][[varietyName]]
              }
            }
          }
        })
      })
      names(resultsLists) <- resultTypes
      resultsLists
    })
    names(resultsByVariety) <- varietyNames
  } else { # leave k out or ordinary, unresampled k-fold cross-validation.
    if(class(results[[1]][["predictions"]]) == "list")
    {
      multipleVarieties <- TRUE
      varietyNames <- names(results[[1]][["predictions"]])
    } else {
      multipleVarieties <- FALSE
      varietyNames <- "None"
    }
    
    resultsByVariety <- lapply(varietyNames, function(varietyName)
    {
      resultsLists <- lapply(resultTypes, function(resultType)
      {
        resultList <- lapply(results, function(sample) sample[[resultType]])
        
        if(resultType == "testSet" || resultType == "predictions" && multipleVarieties == FALSE && class(resultList[[1]]) != "data.frame")
          resultList <- unlist(resultList)
        else if (resultType == "predictions" && multipleVarieties == FALSE) # Predictions are in a data.frame.
          resultList <- do.call(rbind, resultList)
        else if(multipleVarieties == TRUE && resultType == "predictions" && class(resultList[[1]]) != "data.frame")
          unlist(lapply(resultList, "[[", varietyName))
        else if(multipleVarieties == TRUE && resultType == "predictions") # Predictions are factor or numeric.
          do.call(rbind, lapply(resultList, "[[", varietyName))
        else if(multipleVarieties == TRUE && resultType %in% c("ranked", "selected", "tune"))
          lapply(resultList, "[[", varietyName)
        else # For multipleVarieties being FALSE, and not the predicted classes.
          resultList
      })
      names(resultsLists) <- resultTypes
      resultsLists
    })
    names(resultsByVariety) <- varietyNames
  }

  predictionTablesByVariety <- lapply(resultsByVariety, function(resultVariety)
  {
    if(validation == "bootstrap")
    {
      lapply(1:length(results), function(resample)
      {
        switch(class(resultVariety[["predictions"]][[resample]]),
               factor = data.frame(sample = resultVariety[["testSet"]][[resample]], label = resultVariety[["predictions"]][[resample]]),
               numeric = data.frame(sample = resultVariety[["testSet"]][[resample]], score = resultVariety[["predictions"]][[resample]]),
               data.frame = data.frame(sample = resultVariety[["testSet"]][[resample]],
                                       label = resultVariety[["predictions"]][[resample]][, sapply(resultVariety[["predictions"]][[resample]], class) == "factor"],
                                       score = resultVariety[["predictions"]][[resample]][, sapply(resultVariety[["predictions"]][[resample]], class) == "numeric"]))
                
      })
    } else { # leave k out or ordinary, unresampled k-fold cross-validation.
      list(switch(class(resultVariety[["predictions"]]),
             factor = data.frame(sample = resultVariety[["testSet"]], label = resultVariety[["predictions"]]),
             numeric = data.frame(sample = resultVariety[["testSet"]], score = resultVariety[["predictions"]]),
             data.frame = data.frame(sample = resultVariety[["testSet"]],
                                     label = resultVariety[["predictions"]][, sapply(resultVariety[["predictions"]], class) == "factor"],
                                     score = resultVariety[["predictions"]][, sapply(resultVariety[["predictions"]], class) == "numeric"])))
    }
  })

  resultsByVariety <- mapply(function(resultVariety, predictionVariety) # Put formatted prediction tables into results list.
  {
    resultVariety[["predictions"]] <- predictionVariety
    resultVariety
  }, resultsByVariety, predictionTablesByVariety, SIMPLIFY = FALSE)
  
  if(validation == "bootstrap")
  {
    if(bootMode == "fold")
      validationInfo <- list("resampleFold", resamples, folds)
    else
      validationInfo <- list("split", resamples, percent)
  } else if(validation == "leave") {
    validationInfo <- list("leave", leave)
  } else {
    validationInfo <- list("fold", folds)
  }

  if(!"tuneParams" %in% names(params[[match("TrainParams", stagesParamClasses)]]@otherParams)) # No tuning specified.
    resultsByVariety <- lapply(resultsByVariety, function(result) {result[["tune"]] <- list(NULL); result;}) # Shorten it.
  resultsByVariety <- lapply(resultsByVariety, function(result) # Check for empty rankings, and shorten the list.
  {
    if(length(unlist(result[["ranked"]])) == 0)
      result[["ranked"]] <- list(NULL)
    result
  })

  classifyResults <- lapply(varietyNames, function(variety)
  {
    ClassifyResult(datasetName, classificationName, selectParams@selectionName, sampleNames(expression), featureNames(expression),
                   resultsByVariety[[variety]][["ranked"]], resultsByVariety[[variety]][["selected"]], resultsByVariety[[variety]][["predictions"]],
                   pData(expression)[, "class"], validationInfo, resultsByVariety[[variety]][["tune"]])
  })

  names(classifyResults) <- varietyNames
  if(multipleVarieties == FALSE) classifyResults <- classifyResults[[1]]
  classifyResults
})
