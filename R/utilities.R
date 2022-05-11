# Operates on an input data frame, to extract the outcome column(s) and return
# a list with the table of covariates in one element and the outcomes in another.
# The outcomes need to be removed from the data table before predictor training!
.splitDataAndOutcomes <- function(measurements, outcomes, restrict = "numeric")
{ # DataFrame's outcomes variable can be character or factor, so it's a bit involved.
  if(is.character(outcomes) && length(outcomes) > 3 && length(outcomes) != nrow(measurements))
    stop("'outcomes' is a character variable but has more than one element. Either provide a\n",
         "       one to three column names or a factor of the same length as the number of samples.")

  ## String specifies the name of a single outcome column, typically a class.
  if(is.character(outcomes) && length(outcomes) == 1)
  {
    outcomesColumn <- match(outcomes, colnames(measurements))
    if(is.na(outcomesColumn))
      stop("Specified column name of outcomes is not present in the data table.")
    outcomes <- measurements[, classColumn]
    measurements <- measurements[, -classColumn]
    # R version 4 and greater no longer automatically casts character columns to factors because stringsAsFactors
    # is FALSE by default, so it is more likely to be character format these days. Handle it.
    if(class(outcomes) != "factor") # Assume there will be no ordinary regression prediction tasks ... for now.
      outcomes <- factor(outcomes)
  }
  
  # survival's Surv constructor has two inputs for the popular right-censored data and
  # three inputs for less-common interval data.
  if(is.character(outcomes) && length(outcomes) %in% 2:3)
  {
    outcomesColumns <- match(outcomes, colnames(measurements))
    if(any(is.na(outcomesColumns)))
      stop("Specified column names of outcomes is not present in the data table.")
    outcomes <- measurements[, outcomesColumns]
    measurements <- measurements[, -outcomesColumns]
  }
  
  if(is(outcomes, "factor") && length(outcomes) > 3 & length(outcomes) < nrow(measurements))
    stop("The length of outcomes is not equal to the number of samples.")
  
  ## A vector of characters was input by the user. Ensure that it is a factor.
  if(is.character(outcomes) & length(outcomes) == nrow(measurements))
    outcomes <- factor(outcomes)
  
  # Outcomes has columns, so it is tabular. It is inferred to represent survival data.
  if(!is.null(ncol(outcomes)))
  {
    # Assume that event status is in the last column (second for two columns, third for three columns)
    if(!is.null(dim(outcomes)) && length(dim(data)) == 2 && ncol(outcomes) %in% 2:3)
    {
      numberEventTypes <- length(unique(outcomes[, ncol(outcomes)]))
      # Could be one or two kinds of events. All events might be uncensored or censored
      # in a rare but not impossible scenario.
      if(numberEventTypes > 2)
        stop("Number of distinct event types in the last column exceeds 2 but must be 1 or 2.")
      
      
      if(ncol(outcomes) == 2) # Typical, right-censored survival data.
      {
        outcomes <- survival::Surv(outcomes[, 1], outcomes[, 2])
      } else { # Three columns. Therefore, counting process data.
        outcomes <- survival::Surv(outcomes[, 1], outcomes[, 2], outcomes[, 3])
      }
    }
  }
  
  if(!is.null(restrict))
  {
    isDesiredClass <- sapply(measurements, function(column) is(column, restrict))
    measurements <- measurements[, isDesiredClass, drop = FALSE]
    if(ncol(measurements) == 0)
      stop(paste("No features are left after restricting to", restrict, "but at least one must be."))
  }
  
  ###
  # Lets just check that measurements has mcols
  ###
  
  if(is(measurements, "DataFrame")){
    if(is.null(mcols(measurements))){
      message(paste("You have", ncol(measurements), "features and", nrow(measurements), "samples and only one data-type."))
      mcols(measurements)$dataset <- "dataset"
      mcols(measurements)$feature <- colnames(measurements)
  }}
  

  list(measurements = measurements, outcomes = outcomes)
}

# Function to convert a MultiAssayExperiment object into a flat DataFrame table, to enable it
# to be used in typical model building functions.
# Returns a list with a covariate table and and outcomes vector/table, or just a covariate table
# in the case the input is a test data set.
.MAEtoWideTable <- function(measurements, targets = NULL, outcomesColumns = NULL, restrict = "numeric")
{
  if(is.null(targets))
    stop("'targets' is not specified but must be.")
  if(is.null(targets))
    stop("'outcomesColumns' is not specified but must be.")    
  if(!all(targets %in% c(names(measurements), "sampleInfo")))
    stop("Some table names in 'targets' are not assay names in 'measurements' or \"sampleInfo\".")
  sampleInfoColumns <- colnames(MultiAssayExperiment::colData(measurements))
  if(!missing(outcomesColumns) & !all(outcomesColumns %in% sampleInfoColumns))
    stop("Not all column names specified by 'outcomesColumns' found in sample information table.")  

  if("sampleInfo" %in% targets)
  {
    targets <- targets[targets != "sampleInfo"]
    sampleInfoColumnsTrain <- sampleInfoColumns
  } else {
    sampleInfoColumnsTrain <- NULL
  }
  
  if(length(targets) > 0)
  {
    measurements <- measurements[, , targets]
  
    # Get all desired measurements tables and sample information columns (other than the columns representing outcomes).
    # These form the independent variables to be used for making predictions with.
    # Variable names will have names like RNA:BRAF for traceability.
    dataTable <- wideFormat(measurements, colDataCols = union(sampleInfoColumnsTrain, outcomesColumns), check.names = FALSE, collapse = ':')
    rownames(dataTable) <- dataTable[, "primary"]
    S4Vectors::mcols(dataTable)[, "sourceName"] <- gsub("colDataCols", "sampleInfo", S4Vectors::mcols(dataTable)[, "sourceName"])
    colnames(S4Vectors::mcols(dataTable))[1] <- "dataset"
  
    # Sample information variable names not included in column metadata of wide table but only as row names of it.
    # Create a combined column named "feature" which has feature names of the assays as well as the sample information.
    S4Vectors::mcols(dataTable)[, "feature"] <- as.character(S4Vectors::mcols(dataTable)[, "rowname"])
    missingIndices <- is.na(S4Vectors::mcols(dataTable)[, "feature"])
    S4Vectors::mcols(dataTable)[missingIndices, "feature"] <- colnames(dataTable)[missingIndices]
    
    # Finally, a column annotation recording variable name and which table it originated from for all of the source tables.
    S4Vectors::mcols(dataTable) <- S4Vectors::mcols(dataTable)[, c("dataset", "feature")]
  } else { # Must have only been sample information data.
    dataTable <- MultiAssayExperiment::colData(measurements)
  }
  if(!is.null(outcomesColumns)) outcomes <- dataTable[, outcomesColumns]
  
  if(!is.null(restrict))
  {
    isDesiredClass <- sapply(dataTable, function(column) is(column, restrict))
    dataTable <- dataTable[, isDesiredClass, drop = FALSE]
    if(ncol(dataTable) == 0)
      stop(paste("No features are left after restricting to", restrict, "but at least one must be."))
  }

  # Only return independent variables in dataTable for making classifications with.
  # "primary" column is auto-generated by sample information table row names and a duplicate.
  dropColumns <- na.omit(match(c("primary", outcomesColumns), colnames(dataTable)))
  if(length(dropColumns) > 0) dataTable <- dataTable[, -dropColumns]
  
  # Training data table and outcomes for training data.
  if(!is.null(outcomesColumns))
      list(dataTable = dataTable, outcomes = outcomes)
  else # Only test data table for test data input.
    dataTable
}

# For classifiers which use one single function for inputting a training and a testing table,
# and work only for the numeric data type, this checks whether the training and testing tables 
# both have the same set of features and there are at least some numeric features to use,
# after they have been filtered by another function which splits the covariates and the outcomes from the input.
.checkVariablesAndSame <- function(trainingMatrix, testingMatrix)
{
  if(ncol(trainingMatrix) == 0) # Filtering of table removed all columns, leaving nothing to classify with.
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else if(ncol(trainingMatrix) != ncol(testingMatrix))
    stop("Training data set and testing data set contain differing numbers of features.")  
}

# Creates two lists of lists. First has training samples, second has test samples for a range
# of different cross-validation schemes.
#' @import utils
.samplesSplits <- function(crossValParams, outcomes)
{
  if(crossValParams@samplesSplits %in% c("k-Fold", "Permute k-Fold"))
  {
    nPermutations <- ifelse(crossValParams@samplesSplits == "k-Fold", 1, crossValParams@permutations)
    nFolds <- crossValParams@folds
    samplesFolds <- lapply(1:nPermutations, function(permutation)
    {
      # Create maximally-balanced folds, so class balance is about the same in all.
      allFolds <- vector(mode = "list", length = nFolds)
      foldsIndexes <- rep(1:nFolds, length.out = length(outcomes))
      
      foldsIndex = 1
      # Dummy encoding for when outcome is not a class.
      if(is(outcomes, "Surv")) outcomes <- factor(rep("a", length(outcomes)))
      for(outcomeName in levels(outcomes))
      {
        # Permute the indexes of samples in the class.
        whichSamples <- sample(which(outcomes == outcomeName))
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
    # Reorganise into two separate lists, no more nesting.
    list(train = unlist(lapply(samplesFolds, '[[', 1), recursive = FALSE),
         test = unlist(lapply(samplesFolds, '[[', 2), recursive = FALSE))
  } else if(crossValParams@samplesSplits == "Permute Percentage Split") {
    # Take the same percentage of samples from each class to be in training set.
    percent <- crossValParams@percentTest
    samplesTrain <- round((100 - percent) / 100 * table(outcomes))
    samplesTest <- round(percent / 100 * table(outcomes))
    samplesLists <- lapply(1:crossValParams@permutations, function(permutation)
    {
      trainSet <- unlist(mapply(function(outcomeName, number)
      {
        sample(which(outcomes == outcomeName), number)
      }, levels(outcomes), samplesTrain))
      testSet <- setdiff(1:length(classes), trainSet)
      list(trainSet, testSet)
    })
    # Reorganise into two lists: training, testing.
    list(train = lapply(samplesLists, "[[", 1), test = lapply(samplesLists, "[[", 2))
  } else if(crossValParams@samplesSplits == "Leave-k-Out") { # leave k out. 
    testSamples <- as.data.frame(utils::combn(length(outcomes), crossValParams@leave))
    trainingSamples <- lapply(testSamples, function(sample) setdiff(1:length(outcomes), sample))
    list(train = as.list(trainingSamples), test = as.list(testSamples))
  }
}

# Creates a two-column table for tracking the permutation, fold number, or subset of each set
# of test samples (i.e. for leave-k-out scheme).
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

# Add extra variables from within runTest functions to function specified by a params object.
.addIntermediates <- function(params)
{
  intermediateName <- params@intermediate
  intermediates <- list(dynGet(intermediateName, inherits = TRUE))
  if(is.null(names(params@intermediate))) names(intermediates) <- intermediateName else names(intermediates) <- names(params@intermediate)
  params@otherParams <- c(params@otherParams, intermediates)
  params
}


# Carries out one iteration of feature selection. Basically, a ranking function is used to rank
# the features in the training set from best to worst and different top sets are used either for
# predicting on the training set (resubstitution) or nested cross-validation of the training set,
# to find the set of top features which give the best (user-specified) performance measure.
.doSelection <- function(measurementsTrain, outcomesTrain, crossValParams, modellingParams, verbose)
{
  tuneParams <- modellingParams@selectParams@tuneParams
  performanceType <- tuneParams[["performanceType"]]
  topNfeatures <- tuneParams[["nFeatures"]]
  tuneMode <- ifelse("tuneMode" %in% names(tuneParams), tuneParams[["tuneMode"]], crossValParams@tuneMode)
  tuneParams <- tuneParams[-match(c("performanceType", "nFeatures"), names(tuneParams))] # Only used as evaluation metric.
  
  # Make selectParams NULL, since we are currently doing selection and it shouldn't call
  # itself infinitely, but save the parameters of the ranking function for calling the ranking
  # function directly using do.call below.
  featureRanking <- modellingParams@selectParams@featureRanking
  otherParams <- modellingParams@selectParams@otherParams
  doSubset <- modellingParams@selectParams@subsetToSelections 
  modellingParams@selectParams <- NULL
  betterValues <- .ClassifyRenvir[["performanceInfoTable"]][.ClassifyRenvir[["performanceInfoTable"]][, "type"] == performanceType, "better"]
  if(is.function(featureRanking)) # Not a list for ensemble selection.
  {
    paramList <- list(measurementsTrain, outcomesTrain, verbose = verbose)
    paramList <- append(paramList, otherParams) # Used directly by a feature ranking function for rankings of features.
    if(length(tuneParams) == 0) tuneParams <- list(None = "none")
    tuneCombosSelect <- expand.grid(tuneParams, stringsAsFactors = FALSE)

    # Generate feature rankings for each one of the tuning parameter combinations.
    rankings <- lapply(1:nrow(tuneCombosSelect), function(rowIndex)
    {
      tuneCombo <- tuneCombosSelect[rowIndex, , drop = FALSE]
      if(tuneCombo != "none") # Add real parameters before function call.
        paramList <- append(paramList, tuneCombo)
      do.call(featureRanking, paramList)
    })
    
    if(featureRanking@generic == "previousSelection") # Actually selection not ranking.
      return(list(NULL, rankings[[1]], NULL))
    
    if(tuneMode == "none") # Actually selection not ranking.
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
        whichTry <- 1:tuneCombosTrain[rowIndex, "topN"]
        if(doSubset)
        {
          if(is.null(S4Vectors::mcols(measurementsTrain)) ) # There are no different data sets.
          {
            topFeatures <- rankingsVariety[whichTry]
            measurementsTrain <- measurementsTrain[, topFeatures, drop = FALSE] # Features in columns
          } else { # Match to relevant variables, considering data set of them.
            topFeatures <- rankingsVariety[whichTry, ]
            topIDs <-  do.call(paste, topFeatures)
            featuresIDs <- do.call(paste, S4Vectors::mcols(measurementsTrain)[, c("dataset", "feature")])
            topColumns <- match(topIDs, featuresIDs)
            measurementsTrain <- measurementsTrain[, topColumns, drop = FALSE]
          }
        } else { # Pass along features to use.
          modellingParams@trainParams@otherParams <- c(modellingParams@trainParams@otherParams, setNames(list(rankingsVariety[whichTry]), names(modellingParams@trainParams@intermediate)))
        }
        if(ncol(tuneCombosTrain) > 1) # There are some parameters for training.
          modellingParams@trainParams@otherParams <- c(modellingParams@trainParams@otherParams, tuneCombosTrain[rowIndex, 2:ncol(tuneCombosTrain), drop = FALSE])
        modellingParams@trainParams@intermediate <- character(0)
        
        # Do either resubstitution classification or nested-CV classification and calculate the resulting performance metric.
        if(crossValParams@tuneMode == "Resubstitution")
        {
          # Specify measurementsTrain and outcomesTrain for testing, too.
          result <- runTest(measurementsTrain, outcomesTrain, measurementsTrain, outcomesTrain,
                            crossValParams = NULL, modellingParams = modellingParams,
                            verbose = verbose, .iteration = "internal")
          
          predictions <- result[["predictions"]]
          # Classifiers will use a column "class" and survival models will use a column "risk".
          if(class(predictions) == "data.frame")
           predictedOutcomes <- predictions[, na.omit(match(c("class", "risk"), colnames(predictions)))]
          else
           predictedOutcomes <- predictions
          calcExternalPerformance(outcomesTrain, predictedOutcomes, performanceType)
        } else {
           result <- runTests(measurementsTrain, outcomesTrain, crossValParams, modellingParams, verbose = verbose)
           result <- calcCVperformance(result, performanceType)
           median(performance(result)[[performanceType]])
         }
       })

        bestOne <- ifelse(betterValues == "lower", which.min(performances)[1], which.max(performances)[1])
        c(bestOne, performances[bestOne])
      })

      tunePick <- ifelse(betterValues == "lower", which.min(bestPerformers[2, ])[1], which.max(bestPerformers[2, ])[1])
      
      if(verbose == 3)
         message("Features selected.")
      
      tuneRow <- tuneCombosTrain[bestPerformers[1, tunePick], , drop  = FALSE]
      if(ncol(tuneRow) > 1) tuneDetails <- tuneRow[, -1, drop = FALSE] else tuneDetails <- NULL
      
      rankingUse <- rankings[[tunePick]]
      if(is.null(S4Vectors::mcols(measurementsTrain)))
        selection <- rankingUse[1:tuneRow[, "topN"]]
      else # A data frame. Subset the rows.
        selection <- rankingUse[1:tuneRow[, "topN"], ]
      
      list(ranked = rankingUse, selected = selection, tune = tuneDetails)
    } else if(is.list(featureRanking)) { # It is a list of functions for ensemble selection.
      featuresLists <- mapply(function(selector, selParams)
      {
        paramList <- list(measurementsTrain, outcomesTrain, trainParams = trainParams,
                          predictParams = predictParams, verbose = verbose)
        paramList <- append(paramList, selParams)
        do.call(selector, paramList)
      }, modellingParams@selectParams@featureRanking, modellingParams@selectParams@featureRanking, SIMPLIFY = FALSE)

      performances <- sapply(topNfeatures, function(topN)
      {
        topIndices <- unlist(lapply(featuresLists, function(features) features[1:topN]))
        topIndicesCounts <- table(topIndices)
        keep <- names(topIndicesCounts)[topIndicesCounts >= modellingParams@selectParams@minPresence]
        measurementsSelected <- measurementsTrain[, keep, drop = FALSE] # Features in columns
        
        if(crossValParams@tuneMode == "Resubstitution")
        {
          result <- runTest(measurementsSelected, classesTrain,
                            training = 1:nrow(measurementsSelected), testing = 1:nrow(measurementsSelected),
                            crossValParams = NULL, modellingParams,
                            verbose = verbose, .iteration = "internal")
          predictions <- result[["predictions"]]
          if(class(predictions) == "data.frame")
            predictedOutcomes <- predictions[, "class"]
          else
            predictedOutcomes <- predictions
          calcExternalPerformance(outcomesTrain, predictedOutcomes, performanceType)
        } else {
          result <- runTests(measurementsSubset, outcomesTrain, crossValParams, modellingParams, verbose = verbose)
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

# Only for transformations that need to be done within cross-validation.
.doTransform <- function(measurementsTrain, measurementsTest, transformParams, verbose)
{
  paramList <- list(measurementsTrain, measurementsTest)
  if(length(transformParams@otherParams) > 0)
    paramList <- c(paramList, transformParams@otherParams)
  paramList <- c(paramList, verbose = verbose)
  do.call(transformParams@transform, paramList)
}

# Code to create a function call to a training function. Might also do training and testing
# within the same function, so test samples are also passed in case they are needed.
.doTrain <- function(measurementsTrain, outcomesTrain, measurementsTest, outcomesTest, modellingParams, verbose)
{
  tuneChosen <- NULL
  if(!is.null(modellingParams@trainParams@tuneParams) && is.null(modellingParams@selectParams@tuneParams))
  {
    tuneCombos <- expand.grid(modellingParams@trainParams@tuneParams, stringsAsFactors = FALSE)
    modellingParams@trainParams@tuneParams <- NULL
    
    performances <- sapply(1:nrow(tuneCombos), function(rowIndex)
    {
      if(crossValParams@tuneMode == "Resubstitution")
      {
        result <- runTest(measurementsTrain, outcomesTrain, measurementsTest, outcomesTest,
                          crossValParams = NULL, modellingParams,
                          verbose = verbose, .iteration = "internal")
        
        predictions <- result[["predictions"]]
        if(class(predictions) == "data.frame")
          predictedOutcomes <- predictions[, "outcome"]
        else
          predictedOutcomes <- predictions
        calcExternalPerformance(outcomesTest, predictedOutcomes, performanceName)
      } else {
        result <- runTests(measurementsTrain, outcomesTrain,
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
    # Don't name these first two variables. Some classifier functions might use classesTrain and others use outcomesTrain.
    paramList <- list(measurementsTrain, outcomesTrain)
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

# Creates a function call to a prediction function.
.doTest <- function(trained, measurementsTest, predictParams, verbose)
{
  if(!is.null(predictParams@predictor))
  {
    paramList <- list(trained, measurementsTest)
    if(length(predictParams@otherParams) > 0) paramList <- c(paramList, predictParams@otherParams)
    paramList <- c(paramList, verbose = verbose)
    prediction <- do.call(predictParams@predictor, paramList)
  } else { prediction <- trained } # Trained is actually the predictions because only one function, not two.
    if(verbose >= 2)
      message("Prediction completed.")    
    prediction
}

# Converts the characteristics of cross-validation into a pretty character string.
.validationText <- function(crossValParams)
{
  switch(crossValParams@samplesSplits,
  `Permute k-Fold` = paste(crossValParams@permutations, "Permutations,", crossValParams@folds, "Folds"),
  `k-Fold` = paste(crossValParams@folds, "-fold cross-validation", sep = ''),
  `Leave-k-Out` = paste("Leave", crossValParams@leave, "Out"),
  `Permute Percentage Split` = paste(crossValParams@permutations, "Permutations,", crossValParams@percentTest, "% Test"),
  independent = "Independent Set")
}

# Used for ROC area under curve calculation.
# PRtable is a data frame with columns FPR, TPR and class.
# distinctClasses is a vector of all of the class names.
.calcArea <- function(PRtable, distinctClasses)
{
  do.call(rbind, lapply(distinctClasses, function(aClass)
  {
    classTable <- subset(PRtable, class == aClass)
    areaSum <- 0
    for(index in 2:nrow(classTable))
    {
      # Some samples had identical predictions but belong to different classes.
      if(classTable[index, "FPR"] != classTable[index - 1, "FPR"] && classTable[index, "TPR"] != classTable[index - 1, "TPR"])
      {
        newArea <- (classTable[index, "FPR"] - classTable[index - 1, "FPR"]) * classTable[index - 1, "TPR"] + # Rectangle part
         0.5 * (classTable[index, "FPR"] - classTable[index - 1, "FPR"]) * (classTable[index, "TPR"] - classTable[index - 1, "TPR"]) # Triangle part on top.
      } else { # Only one sample with predicted score. Line went either up or right, but not both.
        newArea <- (classTable[index, "FPR"] - classTable[index - 1, "FPR"]) * classTable[index, "TPR"]
      }
      areaSum <- areaSum + newArea
    }
    data.frame(classTable, AUC = round(areaSum, 2), check.names = FALSE)
  }))
}

# Converts features into strings to be displayed in plots.
.getFeaturesStrings <- function(importantFeatures)
{
  # Do nothing if only a simple vector of feature IDs.
  if(!is.null(ncol(importantFeatures[[1]]))) # Data set and feature ID columns.
    importantFeatures <- lapply(importantFeatures, function(features) paste(features[, 1], features[, 2]))
  else if("Pairs" %in% class(importantFeatures[[1]]))
    importantFeatures <- lapply(importantFeatures, function(features) paste(first(features), second(features), sep = '-'))
  importantFeatures
}

# Function to overwrite characteristics which are automatically derived from function names
# by user-specified values.
.filterCharacteristics <- function(characteristics, autoCharacteristics)
{
  # Remove duplication of values for classifiers that have one function for training and 
  # one function for prediction.
  if("Classifier Name" %in% autoCharacteristics[, "characteristic"] && "Predictor Name" %in% autoCharacteristics[, "characteristic"])
  {
    classRow <- which(autoCharacteristics[, "characteristic"] == "Classifier Name")
    predRow <- which(autoCharacteristics[, "characteristic"] == "Predictor Name")
    if(autoCharacteristics[classRow, "value"] == autoCharacteristics[predRow, "value"])
      autoCharacteristics <- autoCharacteristics[-predRow, ]
  }
  # Overwrite automatically-chosen names with user's names.
  if(nrow(autoCharacteristics) > 0 && nrow(characteristics) > 0)
  {
    overwrite <- na.omit(match(characteristics[, "characteristic"], autoCharacteristics[, "characteristic"]))
    if(length(overwrite) > 0)
      autoCharacteristics <- autoCharacteristics[-overwrite, ]
  }
  
  # Merge characteristics tables and return.
  rbind(characteristics, autoCharacteristics)
}

# Don't just plot groups alphabetically, but do so in a meaningful order.
.addUserLevels <- function(plotData, orderingList)
{
  for(orderingIndex in 1:length(orderingList))
  {
    plotData[, names(orderingList)[orderingIndex]] <- factor(plotData[, names(orderingList)[orderingIndex]],
                                                             levels = orderingList[[orderingIndex]])
  }
  plotData
}

# Summary of the features used and the total number of them, no matter if they are a simple type
# or something more complex like Pairs or feature sets.
.summaryFeatures <- function(measurements)
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
  consideredFeatures <- ncol(measurements)
  
  list(allFeatures, featureNames, consideredFeatures)
}

# Function to identify the parameters of an S4 method.
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

# Find the x-axis positions where a set of density functions cross-over.
# The trivial cross-overs at the beginning and end of the data range are removed.
# Used by the mixtures of normals and naive Bayes classifiers.
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

# Samples in the training set are upsampled or downsampled so that the class imbalance is
# removed.
.rebalanceTrainingClasses <- function(measurementsTrain, classesTrain, balancing)
{
  samplesPerClassTrain <- table(classesTrain)
  downsampleTo <- min(samplesPerClassTrain)
  upsampleTo <- max(samplesPerClassTrain)
  trainBalanced <- unlist(mapply(function(classSize, className)
  {
    if(balancing == "downsample" && classSize > downsampleTo)
      sample(which(classesTrain == className), downsampleTo)
    else if(balancing == "upsample" && classSize < upsampleTo)
      sample(which(classesTrain == className), upsampleTo, replace = TRUE)
    else
      which(classesTrain == className)
  }, samplesPerClassTrain, names(samplesPerClassTrain), SIMPLIFY = FALSE))
  measurementsTrain <- measurementsTrain[trainBalanced, ]
  classesTrain <- classesTrain[trainBalanced]
  
  list(measurementsTrain = measurementsTrain, classesTrain = classesTrain)
}

.dlda <- function(x, y, prior = NULL){ # Remove this once sparsediscrim is reinstated to CRAN.
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
    x <- matrix(x, nrow = 1)
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
    posterior <- matrix(posterior, nrow = 1) # Ensure it's always matrix, like just below.
    colnames(posterior) <- names(priors)
  } else {
    posterior <- posterior / rowSums(posterior)
  }
  posterior
}

.dmvnorm_diag <- function(x, mean, sigma) { # Remove once sparsediscrim is reinstated to CRAN.
  exp(sum(dnorm(x, mean=mean, sd=sqrt(sigma), log=TRUE)))
}