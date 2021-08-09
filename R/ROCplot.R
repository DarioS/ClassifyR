setGeneric("ROCplot", function(results, ...)
standardGeneric("ROCplot"))

setMethod("ROCplot", "list", 
          function(results, nBins = sapply(results, totalPredictions),
                   comparisonVariable = c("classificationName", "datasetName", "selectionName", "validation", "None"), lineColours = NULL,
                   lineWidth = 1, fontSizes = c(24, 16, 12, 12, 12), labelPositions = seq(0.0, 1.0, 0.2),
                   plotTitle = "ROC", legendTitle = NULL, xLabel = "False Positive Rate", yLabel = "True Positive Rate",
                   plot = TRUE, showAUC = TRUE)
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")
  if(!requireNamespace("scales", quietly = TRUE))
    stop("The package 'scales' could not be found. Please install it.")   
                      
  comparisonVariable <- match.arg(comparisonVariable)
  ggplot2::theme_set(ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)))
  distinctClasses <- levels(actualClasses(results[[1]]))
  numberDistinctClasses <- length(distinctClasses)

  plotData <- mapply(function(result, resultBins)
  {
    predictions <- result@predictions
    if(class(predictions) == "list") # A list of data.frames. Concatenate them.
      predictions <- do.call(rbind, predictions)

    actualClasses <- actualClasses(result)[match(predictions[, "sample"], sampleNames(result))]
    do.call(rbind, lapply(levels(actualClasses), function(class)
    {
      samplesBins <- .binValues(predictions[, class], resultBins)
      boundaries <- sapply(split(predictions[, class], samplesBins), min)
      totalPositives <- sum(actualClasses == class)
      totalNegatives <- nrow(predictions) - totalPositives
      
      rates <- do.call(rbind, lapply(resultBins:1, function(lowerBin)
      {
        consideredSamples <- samplesBins %in% lowerBin:resultBins
        truePositives <- sum(predictions[consideredSamples, class] >= boundaries[lowerBin] & actualClasses[consideredSamples] == class)
        falsePositives <- sum(predictions[consideredSamples, class] >= boundaries[lowerBin] & actualClasses[consideredSamples] != class)
        TPR <- truePositives / totalPositives
        FPR <- falsePositives / totalNegatives
        data.frame(FPR = FPR, TPR = TPR, class = class)
      }))
      
      rates <- rbind(data.frame(FPR = 0, TPR = 0, class = class), rates)
      
      areaSum = 0
      for(index in 2:nrow(rates))
      {
        areaSum <- areaSum + (rates[index, "FPR"] - rates[index - 1, "FPR"]) * mean(rates[(index - 1) : index, "TPR"])
      }
      validationText <- .validationText(result)                    
      data.frame(dataset = rep(result@datasetName, resultBins + 1),
                 analysis = rep(result@classificationName, resultBins + 1),
                 selection = rep(result@selectResult@selectionName, resultBins + 1),
                 validation = rep(validationText, resultBins + 1),
                 rates,
                 AUC = rep(round(areaSum, 2), resultBins + 1))
    }))
  }, results, nBins, SIMPLIFY = FALSE)
  
  plotData <- do.call(rbind, plotData)
  plotData[, "validation"] <- factor(plotData[, "validation"], levels = unique(plotData[, "validation"])) # Order in which validations were provided.
  
  if(numberDistinctClasses > 2)
    lineColourVariable <- "class"
  else
    lineColourVariable <- comparisonVariable

  if(comparisonVariable != "None")
    if(is.null(lineColours))
      lineColours <- scales::hue_pal()(switch(lineColourVariable, validation = length(unique(plotData[, "validation"])), datasetName = length(unique(plotData[, "dataset"])), classificationName = length(unique(plotData[, "analysis"])), selectionName = length(unique(plotData[, "selection"])), class = numberDistinctClasses))
  if(is.null(legendTitle))
    legendTitle <- switch(lineColourVariable, validation = "Validation", datasetName = "Dataset", classificationName = "Analysis", selectionName = "Selection\nMethod", class = "Class", None = NULL)
  
  if(showAUC == TRUE)
  {
    if(numberDistinctClasses == 2)
    {
      plotData <- subset(plotData, class == distinctClasses[2])
      if(comparisonVariable == "validation")
      {
        plotData[, "validation"] <- paste(plotData[, "validation"], " (AUC", plotData[, "AUC"], ')', sep = '')
        plotData[, "validation"] <- factor(plotData[, "validation"], levels = unique(plotData[, "validation"]))
      } else if(comparisonVariable == "datasetName")
      {
        plotData[, "dataset"] <- paste(plotData[, "dataset"], " (AUC", plotData[, "AUC"], ')', sep = '')
        plotData[, "dataset"] <- factor(plotData[, "dataset"], levels = unique(plotData[, "dataset"]))
      } else if(comparisonVariable == "classificationName")
      {
        plotData[, "analysis"] <- paste(plotData[, "analysis"], " (AUC ", plotData[, "AUC"], ')', sep = '')
        plotData[, "analysis"] <- factor(plotData[, "analysis"], levels = unique(plotData[, "analysis"]))
      } else if(comparisonVariable == "selectionName")
      {
        plotData[, "selection"] <- paste(plotData[, "selection"], " (AUC ", plotData[, "AUC"], ')', sep = '')
        plotData[, "selection"] <- factor(plotData[, "selection"], levels = unique(plotData[, "selection"]))
      }
    } else {
      plotData[, "class"] <- paste(plotData[, "class"], " (AUC ", plotData[, "AUC"], ')', sep = '')
      plotData[, "class"] <- factor(plotData[, "class"], levels = unique(plotData[, "class"]))
    }
  }
  
  plotDataSets <- list(plotData)
  if(numberDistinctClasses > 2)
  {
    comparisonColumn <- switch(comparisonVariable, validation = "validation", selectionName = "selection", datasetName = "dataset", classificationName = "analysis")
    plotDataSets <- split(plotDataSets[[1]], plotDataSets[[1]][, comparisonColumn])
  }

  ROCplots <- lapply(plotDataSets, function(plotData)
              {
                ROCplot <- ggplot2::ggplot(data.frame(plotData), ggplot2::aes_string(x = "FPR", y = "TPR", colour = switch(lineColourVariable, validation = "validation", selectionName = "selection", datasetName = "dataset", classificationName = "analysis", class = "class", None = NULL)), environment = environment()) +
                           ggplot2::geom_line(size = lineWidth) + ggplot2::xlab(NULL) + ggplot2::ylab(NULL) + ggplot2::labs(colour = legendTitle) + ggplot2::geom_segment(x = 0, y = 0, xend = 1, yend = 1, size = lineWidth, colour = "black") + ggplot2::scale_x_continuous(breaks = labelPositions, limits = c(0, 1)) +  ggplot2::scale_y_continuous(breaks = labelPositions, limits = c(0, 1)) +
                           ggplot2::theme(axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), legend.position = c(1, 0), legend.justification = c(1, 0), legend.background = ggplot2::element_rect(fill = "transparent"), legend.title = ggplot2::element_text(size = fontSizes[4], hjust = 0), legend.text = ggplot2::element_text(size = fontSizes[5])) + ggplot2::guides(colour = ggplot2::guide_legend(title.hjust = 0.5)) + ggplot2::scale_colour_manual(values = lineColours)
                
                if(numberDistinctClasses == 2)
                  ROCplot <- ROCplot + ggplot2::xlab(xLabel) + ggplot2::ylab(yLabel) + ggplot2::ggtitle(plotTitle) + ggplot2::theme(axis.title = ggplot2::element_text(size = fontSizes[2]), plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5))
                
                if(length(results) == 1 && showAUC == TRUE && numberDistinctClasses == 2)
                  ROCplot <- ROCplot + ggplot2::annotate("text", x = Inf, y = 0, label = paste("AUC =", plotData[1, "AUC"]), hjust = 1.1, size = fontSizes[2] * 5/14) + ggplot2::theme(legend.position = "none")
                
                if(length(results) > 1 && numberDistinctClasses > 2)
                  ROCplot <- ROCplot + ggplot2::facet_wrap(as.formula(paste('~', comparisonColumn)))
                
                ROCplot
              })

  if(length(ROCplots) == 1)
  {
    ROCplot <- ROCplots[[1]]
  } else {
    ROCplot <- do.call(gridExtra::arrangeGrob,
                       c(ROCplots, nrow = 1, list(top = grid::textGrob(plotTitle, gp = grid::gpar(fontsize = fontSizes[1]), vjust = 0.6),
                                                  left = grid::textGrob(yLabel, gp = grid::gpar(fontsize = fontSizes[2]), rot = 90),
                                                  bottom = grid::textGrob(xLabel, gp = grid::gpar(fontsize = fontSizes[2])))))
  }

  if(plot == TRUE)
    grid::grid.draw(ROCplot)
  
  ROCplot
})
