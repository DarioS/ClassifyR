setGeneric("ROCplot", function(results, ...)
{standardGeneric("ROCplot")})

setMethod("ROCplot", "list", 
          function(results, nBins = 5,
                   lineColourVariable = c("classificationName", "datasetName", "validation", "None"), lineColours = NULL,
                   lineWidth = 1, fontSizes = c(24, 16, 12, 12, 12), labelPositions = seq(0.0, 1.0, 0.2),
                   title = "ROC", xLabel = "False Positive Rate", yLabel = "True Positive Rate",
                   plot = TRUE)
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")
  if(!requireNamespace("scales", quietly = TRUE))
    stop("The package 'scales' could not be found. Please install it.")   
  
  lineColourVariable <- match.arg(lineColourVariable)

  plotData <- lapply(results, function(result)
  {
    predictions <- result@predictions
    if(class(predictions) == "list") # A list of data.frames. Concatenate them.
      predictions <- do.call(rbind, predictions)
    
    actualClasses <- result@actualClasses[predictions[, "sample"]]
    boundaries <- c(min(predictions[, "score"]) - 1, quantile(predictions[, "score"], seq(1/nBins, 1, 1/nBins)))
    samplesBins <- cut(predictions[, "score"], boundaries)
    totalPositives <- sum(actualClasses == levels(actualClasses)[2])
    totalNegatives <- sum(actualClasses == levels(actualClasses)[1])
    rates <- t(sapply(nBins:1, function(lowerBin)
    {
      consideredSamples <- samplesBins %in% levels(samplesBins)[lowerBin:nBins]
      threshold <- levels(samplesBins)[lowerBin:nBins][[1]]
      truePositives <- sum(predictions[consideredSamples, "score"] > boundaries[lowerBin] & actualClasses[consideredSamples] == levels(actualClasses)[2])
      falsePositives <- sum(predictions[consideredSamples, "score"] > boundaries[lowerBin] & actualClasses[consideredSamples] == levels(actualClasses)[1])
      TPR <- truePositives / totalPositives
      FPR <- falsePositives / totalNegatives
      c(FPR = FPR, TPR = TPR)
    }))
    rates <- rbind(rates, data.frame(FPR = 0:1, TPR = 0:1))
    validationText <- .validationText(result)
                        
    data.frame(dataset = rep(result@datasetName, nBins + 2),
               analysis = rep(result@classificationName, nBins + 2),
               validation = rep(validationText, nBins + 2),
               rates)
  })
  plotData <- do.call(rbind, plotData)
  plotData[, "validation"] <- factor(plotData[, "validation"], levels = unique(plotData[, "validation"])) # Order in which validations were provided.
  
  if(lineColourVariable != "None")
    if(is.null(lineColours))
      lineColours <- scales::hue_pal()(switch(lineColourVariable, validation = length(unique(plotData[, "validation"])), datasetName = length(unique(plotData[, "dataset"])), classificationName = length(unique(plotData[, "analysis"]))))

  ROCplot <- ggplot2::ggplot(data.frame(plotData), ggplot2::aes(x = FPR, y = TPR, colour = switch(lineColourVariable, validation = validation, datasetName = dataset, classificationName = analysis, None = NULL)), environment = environment()) +
    ggplot2::geom_line(size = lineWidth) + ggplot2::labs(colour = switch(lineColourVariable, validation = "Validation", datasetName = "Dataset", classificationName = "Analysis", None = NULL)) + ggplot2::geom_segment(x = 0, y = 0, xend = 1, yend = 1, size = lineWidth, colour = "black") + ggplot2::scale_x_continuous(breaks = labelPositions, limits = c(0, 1)) +  ggplot2::scale_y_continuous(breaks = labelPositions, limits = c(0, 1)) + ggplot2::xlab("Top Features") + ggplot2::xlab(xLabel) + ggplot2::ylab(yLabel) +
    ggplot2::ggtitle(title) + ggplot2::theme(axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), legend.position = c(1, 0), legend.justification = c(1, 0), legend.title = ggplot2::element_text(size = fontSizes[4], hjust = 0), legend.text = ggplot2::element_text(size = fontSizes[5]), plot.title = ggplot2::element_text(size = fontSizes[1])) + ggplot2::guides(colour = ggplot2::guide_legend(title.hjust = 0.5))
  if(!is.null(lineColours)) ROCplot <- ROCplot + ggplot2::scale_colour_manual(values = lineColours)
    
  if(plot == TRUE)
    print(ROCplot)
  
  ROCplot
})
