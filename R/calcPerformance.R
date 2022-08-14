#' Add Performance Calculations to a ClassifyResult Object or Calculate for a
#' Pair of Factor Vectors
#' 
#' If \code{calcExternalPerformance} is used, such as when having a vector of
#' known classes and a vector of predicted classes determined outside of the
#' ClassifyR package, a single metric value is calculated. If
#' \code{calcCVperformance} is used, annotates the results of calling
#' \code{\link{runTests}} with one of the user-specified performance measures.
#' 
#' All metrics except Matthews Correlation Coefficient are suitable for
#' evaluating classification scenarios with more than two classes and are
#' reimplementations of those available from Intel DAAL.
#' 
#' If \code{\link{runTests}} was run in resampling mode, one performance
#' measure is produced for every resampling. If the leave-k-out mode was used,
#' then the predictions are concatenated, and one performance measure is
#' calculated for all classifications.
#' 
#' \code{"Balanced Error"} calculates the balanced error rate and is better
#' suited to class-imbalanced data sets than the ordinary error rate specified
#' by \code{"Error"}. \code{"Sample Error"} calculates the error rate of each
#' sample individually. This may help to identify which samples are
#' contributing the most to the overall error rate and check them for
#' confounding factors. Precision, recall and F1 score have micro and macro
#' summary versions. The macro versions are preferable because the metric will
#' not have a good score if there is substantial class imbalance and the
#' classifier predicts all samples as belonging to the majority class.
#' 
#' @aliases calcPerformance calcExternalPerformance calcCVperformance
#' calcExternalPerformance,factor,factor-method calcExternalPerformance,Surv,numeric-method
#' calcCVperformance,ClassifyResult-method
#' @param result An object of class \code{\link{ClassifyResult}}.
#' @param performanceType A character vector of length 1. Default: \code{"Balanced Accuracy"}.
#' Must be one of the following options:
#' \itemize{
#' \item{\code{"Error"}: Ordinary error rate.}
#' \item{\code{"Accuracy"}: Ordinary accuracy.}
#' \item{\code{"Balanced Error"}: Balanced error rate.}
#' \item{\code{"Balanced Accuracy"}: Balanced accuracy.}
#' \item{\code{"Sample Error"}: Error rate for each sample in the data set.}
#' \item{\code{"Sample Accuracy"}: Accuracy for each sample in the data set.}
#' \item{\code{"Micro Precision"}: Sum of the number of correct predictions in
#'         each class, divided by the sum of number of samples in each class.}
#' \item{\code{"Micro Recall"}: Sum of the number of correct predictions in each 
#'         class, divided by the sum of number of samples predicted as
#'         belonging to each class.}
#' \item{\code{"Micro F1"}: F1 score obtained by calculating the
#' harmonic mean of micro precision and micro recall.}
#' \item{\code{"Macro Precision"}: Sum of the ratios of the number of correct predictions
#' in each class to the number of samples in each class, divided by the number of classes.}
#' \item{\code{"Macro Recall"}: Sum of the ratios of the number of correct predictions in each
#' class to the number of samples predicted to be in each class, divided by the number of classes.}
#' \item{\code{"Macro F1"}: F1 score obtained by calculating the harmonic mean of macro precision
#' and macro recall.}
#' \item{\code{"Matthews Correlation Coefficient"}: Matthews Correlation Coefficient (MCC). A score
#' between -1 and 1 indicating how concordant the predicted classes are to the actual classes. Only defined if
#' there are two classes.}
#' \item{\code{"AUC"}: Area Under the Curve. An area ranging from 0 to 1, under the ROC.}
#' \item{\code{"C-index"}: For survival data, the concordance index, for models which produce risk scores. Ranges from 0 to 1.}
#' \item{\code{"Sample C-index"}: Per-individual C-index.}
#' }
#' 
#' @param actualOutcome A factor vector or survival information specifying each sample's known outcome.
#' @param predictedOutcome A factor vector or survival information of the same length as \code{actualOutcome} specifying each sample's predicted outcome.
#' 
#' @return If \code{calcCVperformance} was run, an updated
#' \code{\linkS4class{ClassifyResult}} object, with new metric values in the
#' \code{performance} slot. If \code{calcExternalPerformance} was run, the
#' performance metric value itself.
#' 
#' @author Dario Strbenac
#' @examples
#' 
#'   predictTable <- DataFrame(sample = paste("A", 1:10, sep = ''),
#'                             class = factor(sample(LETTERS[1:2], 50, replace = TRUE)))
#'   actual <- factor(sample(LETTERS[1:2], 10, replace = TRUE))                             
#'   result <- ClassifyResult(DataFrame(characteristic = "Data Set", value = "Example"),
#'                            paste("A", 1:10, sep = ''), list(paste("Gene", 1:50), paste("Gene", 1:50)), list(paste("Gene", 1:5), paste("Gene", 1:10)),
#'                            list(function(oracle){}), NULL, predictTable, actual)
#'   result <- calcCVperformance(result) 
#'   performance(result)
#' 
#' @include classes.R
#' @rdname calcPerformance
#' @usage NULL
#' @export
setGeneric("calcExternalPerformance", function(actualOutcome, predictedOutcome, ...)
standardGeneric("calcExternalPerformance"))

#' @rdname calcPerformance
#' @exportMethod calcExternalPerformance
setMethod("calcExternalPerformance", c("factor", "factor"),
          function(actualOutcome, predictedOutcome, # Both are classes.
                   performanceType = c("Balanced Accuracy", "Balanced Error", "Error", "Accuracy",
                                       "Sample Error", "Sample Accuracy",
                                       "Micro Precision", "Micro Recall",
                                       "Micro F1", "Macro Precision",
                                       "Macro Recall", "Macro F1", "Matthews Correlation Coefficient"))
{
  performanceType <- match.arg(performanceType)
  if(length(levels(actualOutcome)) > 2 && performanceType == "Matthews Correlation Coefficient")
    stop("Error: Matthews Correlation Coefficient specified but data set has more than 2 classes.")
  if(is(predictedOutcome, "factor")) levels(predictedOutcome) <- levels(actualOutcome)
  .calcPerformance(list(actualOutcome), list(predictedOutcome), performanceType = performanceType)[["values"]]
})

#' @rdname calcPerformance
#' @exportMethod calcExternalPerformance
setMethod("calcExternalPerformance", c("Surv", "numeric"),
          function(actualOutcome, predictedOutcome, performanceType = "C-index")
          {
            performanceType <- match.arg(performanceType)
            .calcPerformance(actualOutcome, predictedOutcome, performanceType = performanceType)[["values"]]
          })

#' @rdname calcPerformance
#' @usage NULL
#' @export
setGeneric("calcCVperformance", function(result, ...)
    standardGeneric("calcCVperformance"))

#' @rdname calcPerformance
#' @exportMethod calcCVperformance
setMethod("calcCVperformance", "ClassifyResult",
          function(result, performanceType = c("Balanced Accuracy", "Balanced Error", "Error", "Accuracy",
                                               "Sample Error", "Sample Accuracy",
                                               "Micro Precision", "Micro Recall",
                                               "Micro F1", "Macro Precision",
                                               "Macro Recall", "Macro F1", "Matthews Correlation Coefficient", "AUC", 
                                               "C-index", "Sample C-index"))
{
  performanceType <- match.arg(performanceType)
  actualOutcome <- actualOutcome(result) # Extract the known outcome of each sample.
  
  ### Group by permutation
  if(!performanceType %in% c("Sample Error", "Sample Accuracy"))
  {
    if("permutation" %in% colnames(result@predictions))
      grouping <- result@predictions[, "permutation"]
    else # A set of folds or all leave-k-out predictions or independent train and test sets.
      grouping <- rep(1, nrow(result@predictions))
  }
  
  ### Performance for survival data
  if(performanceType %in% c("C-index", "Sample C-index")) {
    samples <- factor(result@predictions[, "sample"], levels = sampleNames(result))
    performance <- .calcPerformance(actualOutcome = actualOutcome[match(result@predictions[, "sample"], sampleNames(result))],
                                    predictedOutcome = result@predictions[, "risk"], 
                                    samples = samples,
                                    performanceType = performanceType, 
                                    grouping = grouping)
    result@performance[[performance[["name"]]]] <- performance[["values"]]
    return(result)
  }
  
  if(performanceType == "AUC") {
    performance <- .calcPerformance(actualOutcome[match(result@predictions[, "sample"], sampleNames(result))],
                                    result@predictions[, levels(actualOutcome)],
                                    performanceType = performanceType, grouping = grouping)
    result@performance[[performance[["name"]]]] <- performance[["values"]]
    return(result)
  }
  
  ### Performance for data with classes
  if(length(levels(actualOutcome)) > 2 && performanceType == "Matthews Correlation Coefficient")
    stop("Error: Matthews Correlation Coefficient specified but data set has more than 2 classes.")

  classLevels <- levels(actualOutcome)
  samples <- factor(result@predictions[, "sample"], levels = sampleNames(result))
  predictedOutcome <- factor(result@predictions[, "class"], levels = classLevels)
  actualOutcome <- factor(actualOutcome[match(result@predictions[, "sample"], sampleNames(result))], levels = classLevels, ordered = TRUE)
  performance <- .calcPerformance(actualOutcome, predictedOutcome, samples, performanceType, grouping)
  result@performance[[performance[["name"]]]] <- performance[["values"]]
  result
})

#' @importFrom survival concordance
.calcPerformance <- function(actualOutcome, predictedOutcome, samples = NA, performanceType, grouping = NULL)
{
  # Make splitting by group safe.
  if(is(predictedOutcome, "DataFrame")) predictedOutcome <- as.data.frame(predictedOutcome, optional = TRUE)    
  if(performanceType %in% c("Sample Error", "Sample Accuracy"))
  {
    sampleMetricValues <- sapply(levels(samples), function(sampleID)
    {
      consider <- which(samples == sampleID)
      if(performanceType == "Sample Error")
        sum(predictedOutcome[consider] != as.character(actualOutcome[consider]))
      else
        sum(predictedOutcome[consider] == as.character(actualOutcome[consider]))
    })
    performanceValues <- as.numeric(sampleMetricValues / table(samples))
    names(performanceValues) <- levels(samples)
    return(list(name = performanceType, values = performanceValues))
  }

  if(!is.null(grouping))
  {
    actualOutcome <- split(actualOutcome, grouping)
    predictedOutcome <- split(predictedOutcome, grouping)
    if(length(samples) > 1 || !is.na(samples))
    {
      allSamples <- levels(samples)
      samples <- split(samples, grouping)
    }
  }
    
  if(performanceType == "Sample C-index")
  {
    performanceValues <- do.call(rbind, mapply(function(iterationSurv, iterationPredictions, iterationSamples)
    {
      do.call(rbind, lapply(allSamples, function(sampleID)
      {
        sampleIndex <- which(iterationSamples == sampleID)
        otherIndices <- setdiff(seq_along(allSamples), sampleIndex)
        concordants <- discordants <- 0
        iterationSurv <- as.matrix(iterationSurv)
        for(compareIndex in otherIndices)
        {
          if(iterationSurv[sampleIndex, "time"] < iterationSurv[compareIndex, "time"] && iterationPredictions[sampleIndex] > iterationPredictions[compareIndex] && iterationSurv[sampleIndex, "status"] == 1)
          { # Reference sample has shorter time, it is not censored, greater risk. Concordant.
            concordants <- concordants + 1
          } else if(iterationSurv[sampleIndex, "time"] > iterationSurv[compareIndex, "time"] && iterationPredictions[sampleIndex] < iterationPredictions[compareIndex] && iterationSurv[compareIndex, "status"] == 1)
          { # Reference sample has longer time, the comparison sample is not censored, lower risk. Concordant.
            concordants <- concordants + 1
          } else if(iterationSurv[sampleIndex, "time"] < iterationSurv[compareIndex, "time"] && iterationPredictions[sampleIndex] < iterationPredictions[compareIndex] && iterationSurv[sampleIndex, "status"] == 1)
          { # Reference sample has shorter time, it is not censored, but lower risk than comparison sample. Discordant.
            discordants <- discordants + 1
          } else if(iterationSurv[sampleIndex, "time"] > iterationSurv[compareIndex, "time"] && iterationPredictions[sampleIndex] > iterationPredictions[compareIndex] && iterationSurv[compareIndex, "status"] == 1)
          { # Reference sample has longer time, the comparison sample is not censored, but higher risk than comparison sample. Discordant.
            discordants <- discordants + 1
          }
        }
        data.frame(sample = sampleID, concordant = concordants, discordant = discordants)
      }))
    }, actualOutcome, predictedOutcome, samples, SIMPLIFY = FALSE))

    sampleValues <- by(performanceValues[, c("concordant", "discordant")], performanceValues[, "sample"], colSums)
    Cindex <- round(sapply(sampleValues, '[', 1) / (sapply(sampleValues, '[', 1) + sapply(sampleValues, '[', 2)), 2)
    names(Cindex) <- names(sampleValues)
    Cindex <- Cindex[!is.nan(Cindex)] # The individual with the smallest censored time will always have no useful inequalities.
    return(list(name = performanceType, values = Cindex))
  }
  
  if(!is(actualOutcome, "list")) actualOutcome <- list(actualOutcome)
  if(!is(predictedOutcome, "list")) predictedOutcome <- list(predictedOutcome)
  

  if(performanceType %in% c("Accuracy", "Error")) {
    performanceValues <- unlist(mapply(function(iterationClasses, iterationPredictions)
    {
      # Columns are predicted classes, rows are actual classes.
      confusionMatrix <- table(iterationClasses, iterationPredictions)
      totalPredictions <- sum(confusionMatrix)
      correctPredictions <- sum(diag(confusionMatrix))
      diag(confusionMatrix) <- 0
      wrongPredictions <- sum(confusionMatrix)
      if(performanceType == "Accuracy")
        correctPredictions / totalPredictions
      else # It is "error".
        wrongPredictions / totalPredictions
    }, actualOutcome, predictedOutcome, SIMPLIFY = FALSE))
  } else if(performanceType %in% c("Balanced Accuracy", "Balanced Error")) {
    performanceValues <- unlist(mapply(function(iterationClasses, iterationPredictions)
    {
      # Columns are predicted classes, rows are actual classes.
      confusionMatrix <- table(iterationClasses, iterationPredictions)
      classSizes <- rowSums(confusionMatrix)
      classErrors <- classSizes - diag(confusionMatrix)
      if(performanceType == "Balanced Accuracy")
        mean(diag(confusionMatrix) / classSizes)
      else
        mean(classErrors / classSizes)
    }, actualOutcome, predictedOutcome, SIMPLIFY = FALSE))
  } else if(performanceType == "AUC") {
    performanceValues <- unlist(mapply(function(iterationClasses, iterationPredictions)
    {
      classesTable <- do.call(rbind, lapply(levels(iterationClasses), function(class)
      {
        totalPositives <- sum(iterationClasses == class)
        totalNegatives <- sum(iterationClasses != class)
        uniquePredictions <- sort(unique(iterationPredictions[, class]), decreasing = TRUE)
        rates <- do.call(rbind, lapply(uniquePredictions, function(uniquePrediction)
        {
          consideredSamples <- iterationPredictions[, class] >= uniquePrediction
          truePositives <- sum(iterationClasses[consideredSamples] == class)
          falsePositives <- sum(iterationClasses[consideredSamples] != class)
          TPR <- truePositives / totalPositives
          FPR <- falsePositives / totalNegatives
          data.frame(FPR = FPR, TPR = TPR, class = class)
        }))
        rates <- rbind(data.frame(FPR = 0, TPR = 0, class = class), rates)
        rates
      }))
      classesAUC <- .calcArea(classesTable, levels(actualOutcome[[1]]))
      mean(classesAUC[!duplicated(classesAUC[, c("class", "AUC")]), "AUC"]) # Average AUC in iteration.
    }, actualOutcome, predictedOutcome, SIMPLIFY = FALSE))
  } else if(performanceType %in% c("C-index")) {
    performanceValues <- unlist(mapply(function(x, y){
      y <- -y
      survival::concordance(x ~ y)$concordance
    }, actualOutcome, predictedOutcome, SIMPLIFY = FALSE))

    } else { # Metrics for which true positives, true negatives, false positives, false negatives must be calculated.
    performanceValues <- unlist(mapply(function(iterationClasses, iterationPredictions)
    {
      confusionMatrix <- table(iterationClasses, iterationPredictions)
      truePositives <- diag(confusionMatrix)
      falsePositives <- colSums(confusionMatrix) - truePositives
      falseNegatives <- rowSums(confusionMatrix) - truePositives
      trueNegatives <- sum(truePositives) - truePositives
      
      if(performanceType == "Average Accuracy")
      {
        return(sum((truePositives + trueNegatives) / (truePositives + trueNegatives + falsePositives + falseNegatives)) / nrow(confusionMatrix))
      }
      if(performanceType %in% c("Micro Precision", "Micro F1"))
      {
        microP <- sum(truePositives) / sum(truePositives + falsePositives)
        if(performanceType == "Micro Precision") return(microP)
      }
      if(performanceType %in% c("Micro Recall", "Micro F1"))
      {
        microR <- sum(truePositives) / sum(truePositives + falseNegatives)
        if(performanceType == "Micro Recall") return(microR)
      }
      if(performanceType == "Micro F1")
      {
        return(2 * microP * microR / (microP + microR))
      }
      if(performanceType %in% c("Macro Precision", "Macro F1"))
      {
        macroP <- sum(truePositives / (truePositives + falsePositives)) / nrow(confusionMatrix)
        if(performanceType == "Macro Precision") return(macroP)
      }
      if(performanceType %in% c("Macro Recall", "Macro F1"))
      {
        macroR <- sum(truePositives / (truePositives + falseNegatives)) / nrow(confusionMatrix)
        if(performanceType == "Macro Recall") return(macroR)
      }
      if(performanceType == "Macro F1")
      {
        return(2 * macroP * macroR / (macroP + macroR))
      }
      if(performanceType == "Matthews Correlation Coefficient")
      {
        return(unname((truePositives[2] * trueNegatives[2] - falsePositives[2] * falseNegatives[2]) / sqrt((truePositives[2] + falsePositives[2]) * (truePositives[2] + falseNegatives[2]) * (trueNegatives[2] + falsePositives[2]) * (trueNegatives[2] + falseNegatives[2]))))
      }
      
    }, actualOutcome, predictedOutcome, SIMPLIFY = FALSE))
  }

  list(name = performanceType, values = performanceValues)
}