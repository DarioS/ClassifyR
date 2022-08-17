#' Get Frequencies of Feature Selection and Sample-wise Classification Errors
#' 
#' There are two modes. For aggregating feature selection results, the function
#' counts the number of times each feature was selected in all
#' cross-validations. For aggregating classification results, the error rate
#' for each sample is calculated. This is useful in identifying outlier samples
#' that are difficult to classify.
#' 
#' 
#' @aliases distribution distribution,ClassifyResult-method
#' @param result An object of class \code{\link{ClassifyResult}}.
#' @param dataType Whether to calculate sample-wise error rate or the number of
#' times a feature was selected.
#' @param plotType Whether to draw a probability density curve or a histogram.
#' @param summaryType Whether to summarise the feature selections as a
#' percentage or count.
#' @param plot Whether to draw a plot of the frequency of selection or error
#' rate.
#' @param xMax Maximum data value to show in plot.
#' @param fontSizes A vector of length 3. The first number is the size of the
#' title.  The second number is the size of the axes titles. The third number
#' is the size of the axes values.
#' @param ... Further parameters, such as \code{colour} and \code{fill}, passed
#' to \code{\link[ggplot2]{geom_histogram}} or
#' \code{\link[ggplot2]{stat_density}}, depending on the value of
#' \code{plotType}.
#' @return If \code{dataType} is "features", a vector as long as the number of
#' features that were chosen at least once containing the number of times the
#' feature was chosen in cross validations or the percentage of times chosen.
#' If \code{dataType} is "samples", a vector as long as the number of samples,
#' containing the cross-validation error rate of the sample. If \code{plot} is
#' \code{TRUE}, then a plot is also made on the current graphics device.
#' @author Dario Strbenac
#' @examples
#' 
#'   #if(require(sparsediscrim))
#'   #{
#'     data(asthma)
#'     result <- crossValidate(measurements, classes, nRepeats = 5)
#'     featureDistribution <- distribution(result, "features", summaryType = "count",
#'                                         plotType = "histogram", binwidth = 1)
#'     print(head(featureDistribution))
#'   #}
#' @usage NULL
#' @export
setGeneric("distribution", function(result, ...)
           standardGeneric("distribution"))

#' @export
setMethod("distribution", "ClassifyResult", 
          function(result, dataType = c("features", "samples"),
                   plotType = c("density", "histogram"), summaryType = c("percentage", "count"),
                   plot = TRUE, xMax = NULL,
                   fontSizes = c(24, 16, 12), ...)
{
  if(plot == TRUE && !requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")
            
  if(plot == TRUE)
    ggplot2::theme_set(ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)))
            
  dataType <- match.arg(dataType)
  plotType <- match.arg(plotType)
  summaryType <- match.arg(summaryType)
  
  # Automatically choose the axes labels and titles, depending on what is being summarised.
  if(plotType == "density")
    yLabel <- "Density"
  else if(summaryType == "percentage")
    yLabel <- "Percentage"
  else # Count of misclassifications or feature selections.
    yLabel <- "Count"
  if(dataType == "features")
  {
    xLabel <- "Number of Cross-validations"
    title <- "Distribution of Feature Selections"
  } else { # Sample-wise error rate.
    xLabel <- "Number of Misclassifications"
    title <- "Distribution of Sample Misclassifications"
  }

  # summaryType must be percentage if samples are analysed.
  # One CV scheme doesn't guarantee that all samples will be predicted equal number of times.
  if(dataType == "samples") 
  {
    errors <- by(allPredictions, allPredictions[, "sample"], function(samplePredicitons)
              {
                sampleClass <- rep(actualOutcome(result)[samplePredicitons[1, 1]], nrow(samplePredicitons))
                confusion <- table(samplePredicitons[, 2], sampleClass)
                (confusion[upper.tri(confusion)] + confusion[lower.tri(confusion)]) /
                (sum(diag(confusion)) + confusion[upper.tri(confusion)] + confusion[lower.tri(confusion)])
              }) # Sample error rate.
    scores <- rep(NA, length(sampleNames(result)))
    scores[as.numeric(names(errors))] <- errors
    names(scores) <- sampleNames(result)
    scores <- round(scores * 100)
  } else { # features
    chosenFeatures <- chosenFeatureNames(result)
    if(is.vector(chosenFeatures[[1]]))
    {
      allFeatures <- unlist(chosenFeatures)
      allFeaturesText <- allFeatures
    } else if(is(chosenFeatures[[1]], "DataFrame")) {
      allFeatures <- do.call(rbind, chosenFeatures)
      allFeaturesText <- paste(allFeatures[, "assay"], allFeatures[, "feature"], sep = ':')
    } else if("Pairs" %in% class(chosenFeatures[[1]])) {
      allFeatures <- do.call(c, unname(chosenFeatures))
      allFeaturesText <- paste(first(allFeatures), second(allFeatures), sep = ', ')
    } else {
      stop("chosenFeatureNames(result) must be a list of vector, Pairs or DataFrame elements.")
    }
    scores <- table(allFeaturesText)
    if(summaryType == "percentage")
      scores <- scores / length(chosenFeatures) * 100
  }
  
  if(is.null(xMax))
  {
    if(dataType == "features")
      xMax <- max(scores)
    else # Samples
      xMax <- 100 # Error Percentages.
  }
  
  plotData <- data.frame(scores = as.numeric(scores))
  if(plot == TRUE)
  {
    if(!missing(...))
      extras <- list(...)
    else
      extras <- list()
    if(plotType == "density")
    {
      plottedGeom <- do.call(ggplot2::stat_density, c(geom = "path", position = "identity", extras))
    } else { # Histogram plot.
      plottedGeom <- do.call(ggplot2::geom_histogram, extras)
    }
    
    print(ggplot2::ggplot(plotData, ggplot2::aes(x = scores)) + plottedGeom + ggplot2::xlim(0, xMax) +
          ggplot2::xlab(xLabel) + ggplot2::labs(x = xLabel, y = yLabel) + ggplot2::ggtitle(title) +
          ggplot2::theme(axis.title = ggplot2::element_text(size = fontSizes[2]),
                         axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]),
                         plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5)))
  }
  
  # Return scores alongside original chosen features format.
  if(is.vector(chosenFeatures[[1]]))
  { # Simply a vector with names.
    scores
  } else { # Not simple vectors and features. They could be Pairs or data frames.
    isPairs <- "Pairs" %in% class(chosenFeatures[[1]])
    if(isPairs) # Make it DataFrame for counting of the occurrences.
      allFeatures <- as(allFeatures, "DataFrame")
    
    summaryTable <- aggregate(list(count = rep(1, nrow(allFeatures))), as.data.frame(allFeatures, optional = TRUE), length)
    
    if(summaryType == "percentage")
    {
      summaryTable[, 3] <- round(summaryTable[, 3] / length(chosenFeatures) * 100)
      colnames(summaryTable)[3] <- "percentage"
    }
    
    if(isPairs) # Recreate a Pairs object but add metadata of times occurrence.
    {
      pairsSummary <- S4Vectors::Pairs(summaryTable[, "first"], summaryTable[, "second"], summaryTable[, 3])
      colnames(mcols(pairsSummary)) <- colnames(summaryTable[, 3])
      return(pairsSummary)
    } else { # A table of assay and feature.
      summaryTable
    }
  }
})