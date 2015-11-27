setGeneric("distribution", function(result, ...)
           {standardGeneric("distribution")})

setMethod("distribution", "ClassifyResult", 
          function(result, dataType = c("features", "samples"),
                   plotType = c("density", "histogram"), summaryType = c("percentage", "count"),
                   plot = TRUE, xMax = NULL, xLabel = "Percentage of Cross-validations",
                   yLabel = "Density", title = "Distribution of Feature Selections",
                   fontSizes = c(24, 16, 12), ...)
{
  if(plot == TRUE && !requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")            
            
  dataType <- match.arg(dataType)
  plotType <- match.arg(plotType)
  summaryType <- match.arg(summaryType)

  if(dataType == "samples")
  {
    allPredictions <- do.call(rbind, predictions(result))
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
    allFeatures <- unlist(features(result))
    scores <- table(allFeatures)
    names(scores) <- featureNames(result)[as.numeric(names(scores))]
  }
  
  if(dataType == "features" && summaryType == "percentage")
  {
    crossValidations <- length(features(result))
    if(result@validation[[1]] == "fold")
      crossValidations <- crossValidations * length(features(result)[[1]])
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
                         plot.title = ggplot2::element_text(size = fontSizes[1])))
  }
  
  scores
})