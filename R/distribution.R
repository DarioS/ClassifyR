#' Get Frequencies of Feature Selection and Sample Errors
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
#' @param xLabel The label for the x-axis of the plot.
#' @param yLabel The label for the y-axis of the plot.
#' @param title An overall title for the plot.
#' @param fontSizes A vector of length 3. The first number is the size of the
#' title.  The second number is the size of the axes titles. The third number
#' is the size of the axes values.
#' @param ... Further parameters, such as \code{colour} and \code{fill}, passed
#' to \code{\link[ggplot2]{geom_histogram}} or
#' \code{\link[ggplot2]{stat_density}}, depending on the value of
#' \code{plotType}.
#' @return If \code{type} is "features", a vector as long as the number of
#' features that were chosen at least once containing the number of times the
#' feature was chosen in cross validations or the percentage of times chosen.
#' If \code{type} is "samples", a vector as long as the number of samples,
#' containing the cross-validation error rate of the sample. If \code{plot} is
#' \code{TRUE}, then a plot is also made on the current graphics device.
#' @author Dario Strbenac
#' @examples
#' 
#'   #if(require(sparsediscrim))
#'   #{
#'     data(asthma)
#'     CVparams <- CrossValParams(permutations = 5)
#'     result <- runTests(measurements, classes, CVparams, ModellingParams())
#'     featureDistribution <- distribution(result, "features", summaryType = "count",
#'                                         plotType = "histogram",
#'                                         xLabel = "Number of Cross-validations", yLabel = "Count",
#'                                         binwidth = 1)
#'     print(head(featureDistribution))
#'   #}
#' @export
setGeneric("distribution", function(result, ...)
           standardGeneric("distribution"))

setMethod("distribution", "ClassifyResult", 
          function(result, dataType = c("features", "samples"),
                   plotType = c("density", "histogram"), summaryType = c("percentage", "count"),
                   plot = TRUE, xMax = NULL, xLabel = "Percentage of Cross-validations",
                   yLabel = "Density", title = "Distribution of Feature Selections",
                   fontSizes = c(24, 16, 12), ...)
{
  if(plot == TRUE && !requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")
            
  if(plot == TRUE)
    ggplot2::theme_set(ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)))
            
  dataType <- match.arg(dataType)
  plotType <- match.arg(plotType)
  summaryType <- match.arg(summaryType)

  if(dataType == "samples")
  {
    errors <- by(allPredictions, allPredictions[, "sample"], function(samplePredicitons)
             {
                sampleClass <- rep(actualClasses(result)[samplePredicitons[1, 1]], nrow(samplePredicitons))
                confusion <- table(samplePredicitons[, 2], sampleClass)
                (confusion[upper.tri(confusion)] + confusion[lower.tri(confusion)]) /
                (sum(diag(confusion)) + confusion[upper.tri(confusion)] + confusion[lower.tri(confusion)])
             }) # Sample error rate.
    scores <- rep(NA, length(sampleNames(result)))
    scores[as.numeric(names(errors))] <- errors
    names(scores) <- sampleNames(result)
  } else { # features
    chosenFeatures <- features(result)
    if(is.vector(chosenFeatures[[1]])) # No longer numeric row indicies, but character feature IDs.
      allFeatures <- unlist(chosenFeatures)
    else if(is.data.frame(chosenFeatures[[1]])) # Data set and feature ID columns.
      allFeatures <- do.call(rbind, chosenFeatures)
    else if("Pairs" %in% class(chosenFeatures[[1]]))
      allFeatures <- as.data.frame(do.call(c, unname(chosenFeatures)))

    if(is.data.frame(allFeatures))
    {
      if(all(colnames(allFeatures) == c("feature", "dataset")))
        allFeatures <- paste(allFeatures[, "feature"], paste('(', allFeatures[, "dataset"], ')', sep = ''))
      else
        allFeatures <- paste(allFeatures[, "first"], allFeatures[, "second"], sep = ', ')
    }
    scores <- table(allFeatures)
  }

  if(dataType == "features" && summaryType == "percentage")
  {
    crossValidations <- length(features(result))
    scores <- scores / crossValidations * 100
  }
  
  if(is.null(xMax))
  {
    if(dataType == "features")
      xMax <- max(scores)
    else # Samples
      xMax <- 1 # Error rates.
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
  
  scores
})