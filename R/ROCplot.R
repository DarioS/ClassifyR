setGeneric("ROCplot", function(results, ...)
{standardGeneric("ROCplot")})

setMethod("ROCplot", "list", 
          function(results, nBins = sapply(results, totalPredictions),
                   lineColourVariable = c("classificationName", "datasetName", "selectionName", "validation", "None"), lineColours = NULL,
                   lineWidth = 1, fontSizes = c(24, 16, 12, 12, 12), labelPositions = seq(0.0, 1.0, 0.2),
                   plotTitle = "ROC", legendTitle = NULL, xLabel = "False Positive Rate", yLabel = "True Positive Rate",
                   plot = TRUE, showAUC = TRUE)
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")
  if(!requireNamespace("scales", quietly = TRUE))
    stop("The package 'scales' could not be found. Please install it.")   
  lineColourVariable <- match.arg(lineColourVariable)

  plotData <- mapply(function(result, resultBins)
  {
    predictions <- result@predictions
    if(class(predictions) == "list") # A list of data.frames. Concatenate them.
      predictions <- do.call(rbind, predictions)
    
    actualClasses <- result@actualClasses[predictions[, "sample"]]
    samplesBins <- .binValues(predictions[, "score"], resultBins)
    boundaries <- sapply(split(predictions[, "score"], samplesBins), min)
    totalPositives <- sum(actualClasses == levels(actualClasses)[2])
    totalNegatives <- sum(actualClasses == levels(actualClasses)[1])
    rates <- t(sapply(resultBins:1, function(lowerBin)
    {
      consideredSamples <- samplesBins %in% lowerBin:resultBins
      truePositives <- sum(predictions[consideredSamples, "score"] >= boundaries[lowerBin] & actualClasses[consideredSamples] == levels(actualClasses)[2])
      falsePositives <- sum(predictions[consideredSamples, "score"] >= boundaries[lowerBin] & actualClasses[consideredSamples] == levels(actualClasses)[1])
      TPR <- truePositives / totalPositives
      FPR <- falsePositives / totalNegatives
      c(FPR = FPR, TPR = TPR)
    }))
    rates <- rbind(data.frame(FPR = 0, TPR = 0), rates)
    validationText <- .validationText(result)
    areaSum = 0
    for(index in 2:nrow(rates))
    {
      areaSum <- areaSum + (rates[index, "FPR"] - rates[index - 1, "FPR"]) * mean(rates[(index - 1) : index, "TPR"])
    }
                        
    data.frame(dataset = rep(result@datasetName, resultBins + 1),
               analysis = rep(result@classificationName, resultBins + 1),
               selection = rep(result@selectResult@selectionName, resultBins + 1),
               validation = rep(validationText, resultBins + 1),
               rates,
               AUC = rep(round(areaSum, 2), resultBins + 1))
    
  }, results, nBins, SIMPLIFY = FALSE)

  plotData <- do.call(rbind, plotData)
  plotData[, "validation"] <- factor(plotData[, "validation"], levels = unique(plotData[, "validation"])) # Order in which validations were provided.
  
  if(lineColourVariable != "None")
    if(is.null(lineColours))
      lineColours <- scales::hue_pal()(switch(lineColourVariable, validation = length(unique(plotData[, "validation"])), datasetName = length(unique(plotData[, "dataset"])), classificationName = length(unique(plotData[, "analysis"])), selectionName = length(unique(plotData[, "selection"]))))
  if(is.null(legendTitle))
    legendTitle <- switch(lineColourVariable, validation = "Validation", datasetName = "Dataset", classificationName = "Analysis", selectionName = "Selection\nMethod", None = NULL)
  
  if(showAUC == TRUE)
  {
    if(lineColourVariable == "validation")
    {
      plotData[, "validation"] <- paste(plotData[, "validation"], " (AUC", plotData[, "AUC"], ')', sep = '')
      plotData[, "validation"] <- factor(plotData[, "validation"], levels = unique(plotData[, "validation"]))
    } else if(lineColourVariable == "datasetName")
    {
      plotData[, "dataset"] <- paste(plotData[, "dataset"], " (AUC", plotData[, "AUC"], ')', sep = '')
      plotData[, "dataset"] <- factor(plotData[, "dataset"], levels = unique(plotData[, "dataset"]))
    } else if(lineColourVariable == "classificationName")
    {
      plotData[, "analysis"] <- paste(plotData[, "analysis"], " (AUC ", plotData[, "AUC"], ')', sep = '')
      plotData[, "analysis"] <- factor(plotData[, "analysis"], levels = unique(plotData[, "analysis"]))
    } else if(lineColourVariable == "selectionName")
    {
      plotData[, "selection"] <- paste(plotData[, "selection"], " (AUC ", plotData[, "AUC"], ')', sep = '')
      plotData[, "selection"] <- factor(plotData[, "selection"], levels = unique(plotData[, "selection"]))
    }      
  }

  ROCplot <- ggplot2::ggplot(data.frame(plotData), ggplot2::aes_string(x = "FPR", y = "TPR", colour = switch(lineColourVariable, validation = "validation", selectionName = "selection", datasetName = "dataset", classificationName = "analysis", None = NULL)), environment = environment()) +
    ggplot2::geom_line(size = lineWidth) + ggplot2::labs(colour = legendTitle) + ggplot2::geom_segment(x = 0, y = 0, xend = 1, yend = 1, size = lineWidth, colour = "black") + ggplot2::scale_x_continuous(breaks = labelPositions, limits = c(0, 1)) +  ggplot2::scale_y_continuous(breaks = labelPositions, limits = c(0, 1)) + ggplot2::xlab(xLabel) + ggplot2::ylab(yLabel) +
    ggplot2::ggtitle(plotTitle) + ggplot2::theme(axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), legend.position = c(1, 0), legend.justification = c(1, 0), legend.title = ggplot2::element_text(size = fontSizes[4], hjust = 0), legend.text = ggplot2::element_text(size = fontSizes[5]), plot.title = ggplot2::element_text(size = fontSizes[1])) + ggplot2::guides(colour = ggplot2::guide_legend(title.hjust = 0.5)) + ggplot2::scale_colour_manual(values = lineColours)
  
  if(length(results) == 1 && showAUC == TRUE)
      ROCplot <- ROCplot + ggplot2::annotate("text", x = Inf, y = 0, label = paste("AUC =", plotData[1, "AUC"]), hjust = 1.1, size = fontSizes[2]) + ggplot2::theme(legend.position = "none")
    
  if(plot == TRUE)
    print(ROCplot)
  
  ROCplot
})
