.splitDataAndClasses <- function(measurements, classes)
{ # DataFrame methods' class variable can be character or factor, so it's a bit involved.
  if(class(classes) == "character" && length(classes) > 1)
    stop("'classes' is a character variable but has more than one element. Either provide a\n",
         "       single column name or a factor of the same length as the number of samples.")
  
  if(class(classes) == "character")
  {
    classColumn <- match(classes, colnames(measurements))
    if(is.na(classColumn))
      stop("Specified column name of classes is not present in the data table.")
    classes <- measurements[, classColumn]
    measurements <- measurements[, -classColumn]
    if(class(classes) != "factor")
      classes <- factor(classes)
  }
  list(measurements = measurements, classes = classes)
}

.MAEtoWideTable <- function(measurements, targets, restrict = "numeric")
{
  if(is.null(targets))
    stop("'targets' is not specified but must be.")  
  if(!all(targets %in% c(names(measurements), "clinical")))
    stop("Some table names in 'targets' are not assay names in 'measurements' or \"clinical\".")  

  if("clinical" %in% targets)
  {
    clinicalColumns <- colnames(MultiAssayExperiment::colData(measurements))
    targets <- targets[-match("clinical", targets)]
  } else if("class" %in% colnames(MultiAssayExperiment::colData(measurements))) {
    clinicalColumns <- "class"
  } else {
    clinicalColumns <- NULL
  }

  if(length(targets) > 0)
  {
    measurements <- measurements[, , targets]
  
    dataTable <- wideFormat(measurements, colDataCols = clinicalColumns, check.names = FALSE, collapse = ':')
    rownames(dataTable) <- dataTable[, "primary"]
    S4Vectors::mcols(dataTable)[, "sourceName"] <- gsub("colDataCols", "clinical", S4Vectors::mcols(dataTable)[, "sourceName"])
    colnames(S4Vectors::mcols(dataTable))[1] <- "dataset"
  
    S4Vectors::mcols(dataTable)[, "feature"] <- as.character(S4Vectors::mcols(dataTable)[, "rowname"])
    missingIndices <- is.na(S4Vectors::mcols(dataTable)[, "feature"])
    S4Vectors::mcols(dataTable)[missingIndices, "feature"] <- colnames(dataTable)[missingIndices]
    S4Vectors::mcols(dataTable) <- S4Vectors::mcols(dataTable)[, c("dataset", "feature")]
    if("class" %in% colnames(dataTable))
      classes <- dataTable[, "class"]
    else
      classes <- NULL
  } else { # Must have only been clinical data.
    dataTable <- MultiAssayExperiment::colData(measurements)
    classes <- dataTable[, "class"]
  }
  
  if(!is.null(restrict))
  {
    if(restrict == "numeric")
    {
      isNumeric <- sapply(dataTable, is.numeric)
      dataTable <- dataTable[, isNumeric, drop = FALSE]
    } else if(restrict == "integer")
    {
      isInteger <- sapply(dataTable, is.integer)
      dataTable <- dataTable[, isInteger, drop = FALSE]
    }
  }

  # Only return independent variables in the table.
  dropColumns <- na.omit(match(c("primary", "class"), colnames(dataTable)))
  if(length(dropColumns) > 0) dataTable <- dataTable[, -dropColumns]
  
  if(!is.null(classes))
    list(dataTable = dataTable, classes = classes)
  else
    dataTable
}

.checkVariablesAndSame <- function(trainingMatrix, testingMatrix)
{
  if(ncol(trainingMatrix) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else if(ncol(trainingMatrix) != ncol(testingMatrix))
    stop("Training data set and testing data set contain differing numbers of features.")  
}

# Creates two lists of lists. First has training samples, second has test samples for a range
# of different cross-validation schemes.
.samplesSplits <- function(crossValParams, classes)
{
  if(crossValParams@samplesSplits %in% c("k-Fold", "Permute k-Fold"))
  {
    nPermutations <- ifelse(crossValParams@samplesSplits == "k-Fold", 1, crossValParams@permutations)
    nFolds <- crossValParams@folds
    samplesFolds <- lapply(1:nPermutations, function(permutation)
    {
      # Create maximally-balanced folds, so class balance is about the same in all.
      allFolds <- vector(mode = "list", length = nFolds)
      foldsIndexes <- rep(1:nFolds, length.out = length(classes))
      
      foldsIndex = 1
      for(className in levels(classes))
      {
        # Permute the indexes of samples in the class.
        whichSamples <- sample(which(classes == className))
        whichFolds <- foldsIndexes[foldsIndex:(foldsIndex + length(whichSamples) - 1)]
        
        # Put each sample into its fold.
        for(sampleIndex in 1:length(whichSamples))
        {
          allFolds[[whichFolds[sampleIndex]]] <- c(allFolds[[whichFolds[sampleIndex]]], whichSamples[sampleIndex])
        }
        # Move the beginning index to the first new index.
        foldsIndex <- foldsIndex + length(whichSamples)
      }
      
      list(train = lapply(1:nFolds, function(index) unlist(allFolds[setdiff(1:nFolds, index)])),
           test = allFolds
           )
    })
    # Reorganise into three separate lists, no more nesting.
    list(train = unlist(lapply(samplesFolds, '[[', 1), recursive = FALSE),
         test = unlist(lapply(samplesFolds, '[[', 2), recursive = FALSE))
  } else if(crossValParams@samplesSplits == "Permute Percentage Split") {
    # Take the same percentage of samples from each class to be in training set.
    percent <- crossValParams@percentTest
    samplesTrain <- round((100 - percent) / 100 * table(classes))
    samplesTest <- round(percent / 100 * table(classes))
    samplesLists <- lapply(1:crossValParams@permutations, function(permutation)
    {
      trainSet <- unlist(mapply(function(className, number)
      {
        sample(which(classes == className), number)
      }, levels(classes), samplesTrain))
      testSet <- setdiff(1:length(classes), trainSet)
      list(trainSet, testSet)
    })
    # Reorganise into two lists: training, testing.
    list(train = lapply(samplesLists, "[[", 1), test = lapply(samplesLists, "[[", 2))
  } else if(crossValParams@samplesSplits == "Leave-k-Out") { # leave k out. 
    testSamples <- as.data.frame(utils::combn(length(classes), crossValParams@leave))
    trainingSamples <- lapply(testSamples, function(sample) setdiff(1:length(classes), sample))
    list(train = as.list(trainingSamples), test = as.list(testSamples))
  }
}

# Creates a two-column table for tracking the permutation, fold number, or subset of each set
# of test samples.
.splitsTestInfo <- function(crossValParams, samplesSplits)
{
  permutationIDs <- NULL
  foldIDs <- NULL
  subsetIDs <- NULL
  if(crossValParams@samplesSplits %in% c("k-Fold", "Permute k-Fold"))
  {
    foldsSamples <- lengths(samplesSplits[[2]][1:crossValParams@folds])
    totalSamples <- sum(foldsSamples) 
    if(crossValParams@samplesSplits == "Permute k-Fold")
      permutationIDs <- rep(1:crossValParams@permutations, each = totalSamples)
    times <- ifelse(is.null(crossValParams@permutations), 1, crossValParams@permutations)
    foldIDs <- rep(rep(1:crossValParams@folds, foldsSamples), times = times)
  } else if(crossValParams@samplesSplits == "Permute Percentage Split") {
    permutationIDs <- rep(1:crossValParams@permutations, each = length(samplesSplits[[2]][[1]]))
  } else { # Leave-k-out
    totalSamples <- length(unique(unlist(samplesSplits[[2]])))
    subsetIDs <- rep(1:choose(totalSamples, crossValParams@leave), each = crossValParams@leave)
  } 
  
  summaryTable <- cbind(permutation = permutationIDs, fold = foldIDs, subset = subsetIDs)
}

# Add extra variables from within runTest functions to function specified by params.
.addIntermediates <- function(params)
{
  intermediateName <- params@intermediate
  intermediates <- dynGet(intermediateName, inherits = TRUE)
  names(intermediates) <- intermediateName
  params@otherParams <- c(params@otherParams, intermediates)
  params
}

.doSelection <- function(measurements, classes, training, crossValParams, modellingParams, verbose)
{
  names(classes) <- rownames(measurements) # In case training specified by sample IDs rather than numeric indices.
  trainClasses <- droplevels(classes[training])
  tuneParams <- modellingParams@selectParams@tuneParams
  performanceType <- tuneParams[["performanceType"]]
  topNfeatures <- tuneParams[["nFeatures"]]
  tuneParams <- tuneParams[-match(c("performanceType", "nFeatures"), names(tuneParams))] # Only used as evaluation metric.
  featureRanking <- modellingParams@selectParams@featureRanking
  otherParams <- modellingParams@selectParams@otherParams
  modellingParams@selectParams <- NULL
  betterValues <- .ClassifyRenvir[["performanceInfoTable"]][.ClassifyRenvir[["performanceInfoTable"]][, "type"] == performanceType, "better"]
  if(is.function(featureRanking)) # Not a list for ensemble selection.
  {
    paramList <- list(measurements[training, , drop = FALSE], trainClasses, verbose = verbose)
    # Needs training and prediction functions for resubstitution or nested CV performance calculation.
    paramList <- append(paramList, otherParams) # Used directly by a feature selection function for rankings of features.
    if(length(tuneParams) == 0) tuneParams <- list(None = "none")
    tuneCombosSelect <- expand.grid(tuneParams, stringsAsFactors = FALSE)
    
    rankings <- lapply(1:nrow(tuneCombosSelect), function(rowIndex)
    {
      tuneCombo <- tuneCombosSelect[rowIndex, , drop = FALSE]
      if(tuneCombo != "none") # Add real parameters before function call.
        paramList <- append(paramList, tuneCombo)
      do.call(featureRanking, paramList)
    })
    
    if(featureRanking@generic == "previousSelection") # Actually selection not ranking.
      return(list(NULL, rankings[[1]], NULL))
    
    tuneParamsTrain <- list(topN = topNfeatures)
    tuneParamsTrain <- append(tuneParamsTrain, modellingParams@trainParams@tuneParams)
    tuneCombosTrain <- expand.grid(tuneParamsTrain, stringsAsFactors = FALSE)  
    modellingParams@trainParams@tuneParams <- NULL
    bestPerformers <- sapply(rankings, function(rankingsVariety)
    {
      # Creates a matrix. Columns are top n features, rows are varieties (one row if None).
      performances <- sapply(1:nrow(tuneCombosTrain), function(rowIndex)
      {
        topIndices <- rankingsVariety[1:tuneCombosTrain[rowIndex, "topN"]]
        measurementsSelected <- measurements[training, topIndices, drop = FALSE] # Features in columns
        if(ncol(tuneCombosTrain) > 1) # There are some parameters for training.
          modellingParams@trainParams@otherParams <- c(modellingParams@trainParams@otherParams, tuneCombosTrain[rowIndex, 2:ncol(tuneCombosTrain), drop = FALSE])
        
        if(crossValParams@tuneMode == "Resubstitution")
        {
          result <- runTest(measurementsSelected, trainClasses,
                            training = 1:nrow(measurementsSelected), testing = 1:nrow(measurementsSelected),
                            crossValParams = NULL, modellingParams,
                            verbose = verbose, .iteration = "internal")
          
          predictions <- result[["predictions"]]
          if(class(predictions) == "data.frame")
           predictedClasses <- predictions[, "class"]
          else
           predictedClasses <- predictions
          calcExternalPerformance(classes[training], predictedClasses, performanceType)
        } else {
           result <- runTests(measurementsSubset, classes[training], crossValParams, modellingParams, verbose = verbose)
           result <- calcCVperformance(result, performanceType)
           median(performance(aResult)[[performanceType]])
         }
       })
      
        bestOne <- ifelse(betterValues == "lower", which.min(performances)[1], which.max(performances)[1])
        c(bestOne, performances[bestOne])
      })

      tunePick <- ifelse(betterValues == "lower", which.min(bestPerformers[2, ])[1], which.max(bestPerformers[2, ])[1])
        
      if(verbose == 3)
         message("Features selected.")
      
      tuneRow <- tuneCombosTrain[bestPerformers[1, tunePick], , drop  = FALSE]
      if(ncol(tuneRow) > 1) tuneDetails <- tuneRow else tuneDetails <- NULL
      
      list(ranked = rankings[[tunePick]],
           selected = rankings[[tunePick]][1:tuneRow[, "topN"]], tune = tuneDetails)
    } else if(is.list(featureRanking)) { # It is a list of functions for ensemble selection.
      featuresLists <- mapply(function(selector, selParams)
      {
        paramList <- list(measurements[training, , drop = FALSE], trainClasses, trainParams = trainParams,
                          predictParams = predictParams, verbose = verbose)
        paramList <- append(paramList, selParams)
        do.call(selector, paramList)
      }, modellingParams@selectParams@featureRanking, modellingParams@selectParams@featureRanking, SIMPLIFY = FALSE)

      performances <- sapply(topNfeatures, function(topN)
      {
        topIndices <- unlist(lapply(featuresLists, function(features) features[1:topN]))
        topIndicesCounts <- table(topIndices)
        keep <- names(topIndicesCounts)[topIndicesCounts >= modellingParams@selectParams@minPresence]
        measurementsSelected <- measurements[training, keep, drop = FALSE] # Features in columns
        
        if(crossValParams@tuneMode == "Resubstitution")
        {
          result <- runTest(measurementsSelected, trainClasses,
                            training = 1:nrow(measurementsSelected), testing = 1:nrow(measurementsSelected),
                            crossValParams = NULL, modellingParams,
                            verbose = verbose, .iteration = "internal")
          predictions <- result[["predictions"]]
          if(class(predictions) == "data.frame")
            predictedClasses <- predictions[, "class"]
          else
            predictedClasses <- predictions
          calcExternalPerformance(classes[training], predictedClasses, performanceType)
        } else {
          result <- runTests(measurementsSubset, classes[training], crossValParams, modellingParams, verbose = verbose)
          result <- calcCVperformance(result, performanceType)
          median(performance(aResult)[[performanceType]])
        }
      })
      bestOne <- ifelse(betterValues == "lower", which.min(performances)[1], which.max(performances)[1])
      
      selectedFeatures <- unlist(lapply(featuresLists, function(featuresList) featuresList[1:topNfeatures[bestOne]]))
      names(table(selectedFeatures))[table(selectedFeatures) >= modellingParams@selectParams@minPresence]
      
      list(NULL, selectedFeatures, NULL)
    } else { # Previous selection
      selectedFeatures <- 
      list(NULL, selectedFeatures, NULL)
    }
}

.doTransform <- function(measurements, training, transformParams, verbose)
{
  paramList <- list(measurements, training = training) # Often a central point, like the mean, is used for subtraction or standardisation of values. Pass this to the transformation function.
  if(length(transformParams@otherParams) > 0)
    paramList <- c(paramList, transformParams@otherParams)
  paramList <- c(paramList, verbose = verbose)
  do.call(transformParams@transform, paramList)
}

.doTrain <- function(measurements, classes, training, testing, modellingParams, verbose)
{
  names(classes) <- rownames(measurements) # In case training or testing specified by sample IDs rather than numeric indices.
  trainClasses <- droplevels(classes[training])
  measurementsTrain <- measurements[training, , drop = FALSE]
  measurementsTest <- measurements[testing, , drop = FALSE]
  
  tuneChosen <- NULL
  if(!is.null(modellingParams@trainParams@tuneParams) && is.null(modellingParams@selectParams@tuneParams))
  {
    tuneCombos <- expand.grid(modellingParams@trainParams@tuneParams, stringsAsFactors = FALSE)
    modellingParams@trainParams@tuneParams <- NULL
    
    performances <- sapply(1:nrow(tuneCombos), function(rowIndex)
    {
      if(crossValParams@tuneMode == "Resubstitution")
      {
        result <- runTest(measurements, classes,
                          training = training, testing = testing,
                          crossValParams = NULL, modellingParams,
                          verbose = verbose, .iteration = "internal")
        
        predictions <- result[["predictions"]]
        if(class(predictions[[1]]) == "data.frame")
          predictedClasses <- lapply(predictions, function(set) set[, "class"])
        else
          predictedClasses <- predictions
        sapply(predictedClasses, function(classSet) calcExternalPerformance(classes, classSet, performanceName))
      } else {
        result <- runTests(measurementsSubset, classes,
                           crossValParams, modellingParams,
                           verbose = verbose, .iteration = "internal")
        result <- calcCVperformance(result, performanceName)
        median(predictions(result)["performanceType"])
      }
    })
    betterValues <- .ClassifyRenvir[["performanceInfoTable"]][.ClassifyRenvir[["performanceInfoTable"]][, "type"] == performanceType, "better"]
    bestOne <- ifelse(betterValues == "lower", which.min(performances)[1], which.max(performances)[1])
    tuneChosen <- tuneCombos[bestOne, , drop = FALSE]
    modellingParams@trainParams@otherParams <- tuneChosen
  }

  if(modellingParams@trainParams@classifier@generic != "previousTrained")
    paramList <- list(measurementsTrain, trainClasses)
  else # Don't pass the measurements and classes, because a pre-existing classifier is used.
    paramList <- list()
  if(is.null(modellingParams@predictParams)) # One function does both training and testing.
      paramList <- c(paramList, measurementsTest)
    
  if(length(modellingParams@trainParams@otherParams) > 0)
    paramList <- c(paramList, modellingParams@trainParams@otherParams)
  paramList <- c(paramList, verbose = verbose)

  trained <- do.call(modellingParams@trainParams@classifier, paramList)
  if(verbose >= 2)
    message("Training completed.")  
  
  list(model = trained, tune = tuneChosen)
}

.doTest <- function(trained, measurements, testing, predictParams, verbose)
{
  if(!is.null(predictParams@predictor))
  {
    testMeasurements <- measurements[testing, , drop = FALSE]
  
    paramList <- list(trained, testMeasurements)
    if(length(predictParams@otherParams) > 0) paramList <- c(paramList, predictParams@otherParams)
      paramList <- c(paramList, verbose = verbose)
      prediction <- do.call(predictParams@predictor, paramList)
    } else { prediction <- trained } # Trained is actually the predictions because only one function, not two.
    
    if(verbose >= 2)
      message("Prediction completed.")    
    prediction
}

.validationText <- function(crossValParams)
{
  switch(crossValParams@samplesSplits,
  `Permute k-Fold` = paste(crossValParams@permutations, "Permutations,", crossValParams@folds, "Folds"),
  `k-Fold` = paste(crossValParams@folds, "-fold cross-validation", sep = ''),
  `Leave-k-Out` = paste("Leave", crossValParams@leave, "Out"),
  `Permute Percentage Split` = paste(crossValParams@permutations, "Permutations,", crossValParams@percentTest, "% Test"),
  independent = "Independent Set")
}

.binValues <- function(values, nBins)
{
  ordering <- order(values)
  binID <- rep(1:nBins, each = length(values) / nBins)
  if(length(binID) < length(values))
    binID <- c(binID, rep(max(binID) + 1, length(values) - length(binID)))
  bins <- split(ordering, binID)  
  binID <- numeric()
  binID[unlist(bins)] <- rep(as.numeric(names(bins)), sapply(bins, length))
  binID
}

.getFeaturesStrings <- function(importantFeatures)
{
  if(is.data.frame(importantFeatures[[1]])) # Data set and feature ID columns.
    importantFeatures <- lapply(importantFeatures, function(features) paste(features[, 1], features[, 2]))
  else if("Pairs" %in% class(importantFeatures[[1]]))
    importantFeatures <- lapply(importantFeatures, function(features) paste(first(features), second(features)))
  else if(class(importantFeatures[[1]]) == "list" && is.character(importantFeatures[[1]][[1]])) # Two-level list, such as generated by permuting and folding.
    importantFeatures <- unlist(importantFeatures, recursive = FALSE)
  else if(class(importantFeatures[[1]]) == "list" && is.data.frame(importantFeatures[[1]][[1]])) # Data set and feature ID columns.
    importantFeatures <- unlist(lapply(importantFeatures, function(folds) lapply(folds, function(fold) paste(fold[, 1], fold[, 2])), recursive = FALSE))
  else if(class(importantFeatures[[1]]) == "list" && "Pairs" %in% class(importantFeatures[[1]]))
    importantFeatures <- unlist(lapply(importantFeatures, function(folds) lapply(folds, function(fold) paste(first(fold), second(fold))), recursive = FALSE))  
  importantFeatures
}

.filterCharacteristics <- function(characteristics, autoCharacteristics)
{
  if("Classifier Name" %in% autoCharacteristics[, "characteristic"] && "Predictor Name" %in% autoCharacteristics[, "characteristic"])
  {
    classRow <- which(autoCharacteristics[, "characteristic"] == "Classifier Name")
    predRow <- which(autoCharacteristics[, "characteristic"] == "Predictor Name")
    if(autoCharacteristics[classRow, "value"] == autoCharacteristics[predRow, "value"])
      autoCharacteristics <- autoCharacteristics[-predRow, ]
  }
  if(nrow(autoCharacteristics) > 0 && nrow(characteristics) > 0)
  {
    overwrite <- na.omit(match(characteristics[, "characteristic"], autoCharacteristics[, "characteristic"]))
    if(length(overwrite) > 0)
      autoCharacteristics <- autoCharacteristics[-overwrite, ]
  }
  characteristics <- rbind(characteristics, autoCharacteristics)
}

.addUserLevels <- function(plotData, orderingList) # Don't just plot groups alphabetically but meaningful order.
{
  for(orderingIndex in 1:length(orderingList))
  {
    plotData[, names(orderingList)[orderingIndex]] <- factor(plotData[, names(orderingList)[orderingIndex]],
                                                             levels = orderingList[[orderingIndex]])
  }
  plotData
}

.summaryFeatures <- function(measurements, prevalidate)
{
  # MultiAssayExperiment has feature details in mcols.
  if(!is.null(S4Vectors::mcols(measurements)))
  {
    allFeatures <- S4Vectors::mcols(measurements)
    featureNames <- S4Vectors::mcols(measurements)[, "feature"]
  } else {
    allFeatures <- colnames(measurements)
    featureNames <- colnames(measurements)
  }
  
  # Could refer to features or feature sets, depending on if a selection method utilising feature sets is used.
  if(prevalidate == FALSE)
  {
    consideredFeatures <- ncol(measurements)
  } else {
    consideredFeatures <- nrow(subset(mcols(measurements), dataset == "clinical")) +
      length(unique(subset(mcols(measurements), dataset != "clinical")[, "dataset"]))
  }
  list(allFeatures, featureNames, consideredFeatures)
}

.methodFormals <- function(f, signature) {
  tryCatch({
    fdef <- getGeneric(f)
    method <- selectMethod(fdef, signature)
    genFormals <- base::formals(fdef)
    b <- body(method)
    if(is(b, "{") && is(b[[2]], "<-") && identical(b[[2]][[2]], as.name(".local"))) {
      local <- eval(b[[2]][[3]])
      if(is.function(local))
        return(formals(local))
      warning("Expected a .local assignment to be a function. Corrupted method?")
    }
    genFormals
  },
    error = function(error) {
      formals(f)
    })
}

.densitiesCrossover <- function(densities) # A list of densities created by splinefun.
{
  
  if(!all(table(unlist(lapply(densities, function(density) density[['x']]))) == length(densities)))
    stop("x positions are not the same for all of the densities.")
  
  lapply(1:length(densities), function(densityIndex) # All crossing points with other class densities.
  {
    unlist(lapply(setdiff(1:length(densities), densityIndex), function(otherIndex)
    {
      allDifferences <- densities[[densityIndex]][['y']] - densities[[otherIndex]][['y']]
      crosses <- which(diff(sign(allDifferences)) != 0)
      crosses <- sapply(crosses, function(cross) # Refine location for plateaus.
      {
        isSmall <- rle(allDifferences[(cross+1):length(allDifferences)] < 0.000001)
        if(isSmall[["values"]][1] == "TRUE")
          cross <- cross + isSmall[["lengths"]][1] / 2
        cross
      })
      if(length(crosses) > 1 && densities[[densityIndex]][['y']][crosses[1]] < 0.000001 && densities[[densityIndex]][['y']][crosses[length(crosses)]] < 0.000001)
        crosses <- crosses[-c(1, length(crosses))] # Remove crossings at ends of densities.      
      densities[[densityIndex]][['x']][crosses]
    }))
  })
}

.dlda <- function(x, y, prior = NULL){ # Remove once sparsediscrim is reinstated to CRAN.
  obj <- list()
  obj$labels <- y
  obj$N <- length(y)
  obj$p <- ncol(x)
  obj$groups <- levels(y)
  obj$num_groups <- nlevels(y)

  est_mean <- "mle"

  # Error Checking
  if (!is.null(prior)) {
    if (length(prior) != obj$num_groups) {
      stop("The number of 'prior' probabilities must match the number of classes in 'y'.")
    }
    if (any(prior <= 0)) {
      stop("The 'prior' probabilities must be nonnegative.")
    }
    if (sum(prior) != 1) {
      stop("The 'prior' probabilities must sum to one.")
    }
  }
  if (any(table(y) < 2)) {
    stop("There must be at least 2 observations in each class.")
  }

  # By default, we estimate the 'a priori' probabilities of class membership with
  # the MLEs (the sample proportions).
  if (is.null(prior)) {
    prior <- as.vector(table(y) / length(y))
  }

  # For each class, we calculate the MLEs (or specified alternative estimators)
  # for each parameter used in the DLDA classifier. The 'est' list contains the
  # estimators for each class.
  obj$est <- tapply(seq_along(y), y, function(i) {
    stats <- list()
    stats$n <- length(i)
    stats$xbar <- colMeans(x[i, , drop = FALSE])
    stats$var <- with(stats, (n - 1) / n * apply(x[i, , drop = FALSE], 2, var))
    stats
  })

  # Calculates the pooled variance across all classes.
  obj$var_pool <- Reduce('+', lapply(obj$est, function(x) x$n * x$var)) / obj$N

  # Add each element in 'prior' to the corresponding obj$est$prior
  for(k in seq_len(obj$num_groups)) {
    obj$est[[k]]$prior <- prior[k]
  }
  class(obj) <- "dlda"
  obj
}

.predict <- function(object, newdata, ...) { # Remove once sparsediscrim is reinstated to CRAN.
  if (!inherits(object, "dlda"))  {
    stop("object not of class 'dlda'")
  }
  if (is.vector(newdata)) {
    newdata <- as.matrix(newdata)
  }

  scores <- apply(newdata, 1, function(obs) {
    sapply(object$est, function(class_est) {
      with(class_est, sum((obs - xbar)^2 / object$var_pool) + log(prior))
    })
  })

  if (is.vector(scores)) {
    min_scores <- which.min(scores)
  } else {
    min_scores <- apply(scores, 2, which.min)
  }

  # Posterior probabilities via Bayes Theorem
  means <- lapply(object$est, "[[", "xbar")
  covs <- replicate(n=object$num_groups, object$var_pool, simplify=FALSE)
  priors <- lapply(object$est, "[[", "prior")
  posterior <- .posterior_probs(x=newdata,
                               means=means,
                               covs=covs,
                               priors=priors)

  class <- factor(object$groups[min_scores], levels = object$groups)

  list(class = class, scores = scores, posterior = posterior)
}

.posterior_probs <- function(x, means, covs, priors) { # Remove once sparsediscrim is reinstated to CRAN.
  if (is.vector(x)) {
    x <- matrix(x, nrow=1)
  }
  x <- as.matrix(x)

  posterior <- mapply(function(xbar_k, cov_k, prior_k) {
    if (is.vector(cov_k)) {
      post_k <- apply(x, 1, function(obs) {
        .dmvnorm_diag(x=obs, mean=xbar_k, sigma=cov_k)
      })
    } else {
      post_k <- dmvnorm(x=x, mean=xbar_k, sigma=cov_k)
    }
    prior_k * post_k
  }, means, covs, priors)

  if (is.vector(posterior)) {
    posterior <- posterior / sum(posterior)
  } else {
    posterior <- posterior / rowSums(posterior)
  }

  posterior
}

.dmvnorm_diag <- function(x, mean, sigma) { # Remove once sparsediscrim is reinstated to CRAN.
  exp(sum(dnorm(x, mean=mean, sd=sqrt(sigma), log=TRUE)))
}

.rebalanceTrainingClasses <- function(measurements, classes, training, testing, balancing)
{
  samplesPerClassTrain <- table(classes[training])
  downsampleTo <- min(samplesPerClassTrain)
  upsampleTo <- max(samplesPerClassTrain)
  trainBalanced <- as.vector(mapply(function(classSize, className)
  {
    if(balancing == "downsample" && classSize > downsampleTo)
      sample(which(classes[training] == className), downsampleTo)
    else if(balancing == "upsample" && classSize < upsampleTo)
      sample(which(classes[training] == className), upsampleTo, replace = TRUE)
    else
      which(classes[training] == className)
  }, samplesPerClassTrain, names(samplesPerClassTrain)))
  measurements <- measurements[c(training[trainBalanced], testing), ]
  classes <- classes[c(training[trainBalanced], testing)]
  training <- 1:length(trainBalanced)
  testing <- (length(training) + 1):length(classes)
  
  list(measurements = measurements, classes = classes, training = training, testing = testing)
}