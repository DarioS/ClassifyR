setGeneric("performancePlot", function(results, ...)
{standardGeneric("performancePlot")})

setMethod("performancePlot", "list", 
          function(results,
                   aggregate = character(),
                   xVariable = c("classificationName", "datasetName", "selectionName", "validation"),
                   performanceName = NULL,
                   boxFillColouring = c("classificationName", "datasetName", "selectionName", "validation", "None"),
                   boxFillColours = NULL,
                   boxLineColouring = c("classificationName", "datasetName", "selectionName", "validation", "None"),
                   boxLineColours = NULL,
                   rowVariable = c("None", "validation", "datasetName", "classificationName", "selectionName"),
                   columnVariable = c("datasetName", "classificationName", "validation", "selectionName", "None"),
                   yMax = 1, fontSizes = c(24, 16, 12, 12), title = NULL,
                   xLabel = "Analysis", yLabel = performanceName,
                   margin = grid::unit(c(0, 1, 1, 0), "lines"), rotate90 = FALSE, showLegend = TRUE, plot = TRUE)
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")             
  if(!requireNamespace("scales", quietly = TRUE))
    stop("The package 'scales' could not be found. Please install it.")
  
  xVariable <- match.arg(xVariable)
  boxFillColouring <- match.arg(boxFillColouring)
  boxLineColouring <- match.arg(boxLineColouring)
  rowVariable <- match.arg(rowVariable)
  columnVariable <- match.arg(columnVariable)  
  if(is.null(performanceName)) stop("Please specify a performance measure to plot.")
  performances <- lapply(results, function(result)
                         {
                           if(!performanceName %in% names(result@performance))
                             stop("Performance measure not calculated for ", result@classificationName)
                           if(xVariable == "selectionName") slotValue <- slot(slot(result, "selectResult"), "selectionName")
                           else slotValue <- slot(result, xVariable)
                           if(slotValue %in% aggregate) mean(result@performance[[performanceName]])
                           else result@performance[[performanceName]]
  })
  performanceLengths <- sapply(performances, length)
  
  plotData <- data.frame(dataset = rep(sapply(results, function(result) result@datasetName), performanceLengths),
                         analysis = rep(sapply(results, function(result) result@classificationName), performanceLengths),
                         selection = rep(sapply(results, function(result) result@selectResult@selectionName), performanceLengths),
                         validation = rep(sapply(results, function(result) .validationText(result)), performanceLengths),
                         performance = unlist(performances))
  
  if(boxFillColouring != "None")
    if(is.null(boxFillColours)) boxFillColours <- scales::hue_pal()(switch(boxFillColouring, validation = length(unique(plotData[, "validation"])), datasetName = length(unique(plotData[, "dataset"])), classificationName = length(unique(plotData[, "analysis"])), selectionName = length(unique(plotData[, "selection"]))))
  if(boxLineColouring != "None")
    if(is.null(boxLineColours)) boxLineColours <- scales::hue_pal(direction = -1)(switch(boxLineColouring, validation = length(unique(plotData[, "validation"])), datasetName = length(unique(plotData[, "dataset"])), classificationName = length(unique(plotData[, "analysis"])), selectionName = length(unique(plotData[, "selection"]))))
  
  plotData[, "dataset"] <- factor(plotData[, "dataset"], levels = unique(plotData[, "dataset"]))
  plotData[, "analysis"] <- factor(plotData[, "analysis"], levels = unique(plotData[, "analysis"]))
  plotData[, "selection"] <- factor(plotData[, "selection"], levels = unique(plotData[, "selection"]))
  plotData[, "validation"] <- factor(plotData[, "validation"], levels = unique(plotData[, "validation"])) # Order in which validations were provided.
  
  if(rotate90 == TRUE)
  {
    switch(xVariable, validation = plotData[, "validation"] <- factor(plotData[, "validation"], levels = rev(levels(plotData[, "validation"]))),
           datasetName = plotData[, "dataset"] <- factor(plotData[, "dataset"], levels = rev(levels(plotData[, "dataset"]))),
           selectionName = plotData[, "selection"] <- factor(plotData[, "selection"], levels = rev(levels(plotData[, "selection"]))),
           classificationName = plotData[, "analysis"] <- factor(plotData[, "analysis"], levels = rev(levels(plotData[, "analysis"]))))
  }

  legendPosition <- ifelse(showLegend == TRUE, "right", "none")
  performancePlot <- ggplot2::ggplot() + ggplot2::scale_y_continuous(limits = c(0, yMax)) + ggplot2::xlab(xLabel) + ggplot2::ylab(yLabel) +
                          ggplot2::ggtitle(title) + ggplot2::theme(legend.position = legendPosition, axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), plot.title = ggplot2::element_text(size = fontSizes[1]), plot.margin = margin)

  # Compatible with ggplot2 1.0.1.
  if(!is.null(boxFillColours))
    performancePlot <- performancePlot + ggplot2::scale_fill_manual(values = boxFillColours)
  if(!is.null(boxLineColours))
    performancePlot <- performancePlot + ggplot2::scale_fill_manual(values = boxLineColours)
  
  performanceCounts <- as.data.frame(table(plotData[, 1:4]))
  if(any(performanceCounts[, "Freq"] > 1))
  {
    multiplePerformance <- subset(performanceCounts, Freq > 1)
    multipleRows <- rowSums(apply(multiplePerformance[, 1:4], 1, function(parameters)
      apply(plotData[, 1:4], 1, function(plotRow) all(plotRow == parameters))
    )) > 0
    multiPlotData <- plotData[multipleRows, ]
    performancePlot <- performancePlot + ggplot2::geom_boxplot(data = multiPlotData, ggplot2::aes_string(x = switch(xVariable, validation = "validation", datasetName = "dataset", classificationName = "analysis", selectionName = "selection"), y = "performance",
                                                                                                         fill = switch(boxFillColouring, validation = "validation", datasetName = "dataset", classificationName = "analysis", selectionName = "selection", None = NULL), colour = switch(boxLineColouring, validation = "validation", datasetName = "dataset", classificationName = "analysis", selectionName = "selection", None = NULL)))
  }
  if(any(performanceCounts[, "Freq"] == 1))
  {
    singlePerformance <- subset(performanceCounts, Freq == 1)
    singleRows <- rowSums(apply(singlePerformance[, 1:4], 1, function(parameters)
      apply(plotData[, 1:4], 1, function(plotRow) all(plotRow == parameters))
    )) > 0
    singlePlotData <- plotData[singleRows, ]
    performancePlot <- performancePlot + ggplot2::geom_bar(data = singlePlotData, stat = "identity", ggplot2::aes_string(x = switch(xVariable, validation = "validation", datasetName = "dataset", classificationName = "analysis", selectionName = "selection"), y = "performance", fill = switch(boxFillColouring, validation = "validation", datasetName = "dataset", classificationName = "analysis", selectionName = "selection", None = NULL), colour = switch(boxLineColouring, validation = "validation", datasetName = "dataset", selectionName = "selection", classificationName = "analysis", None = NULL)))    
  }
  
  if(rotate90 == TRUE)
    performancePlot <- performancePlot + ggplot2::coord_flip()
  if(legendPosition != "none")
    performancePlot <- performancePlot + ggplot2::labs(colour = switch(boxLineColouring, validation = "Validation", datasetName = "Dataset", classificationName = "Analysis", classificationName = "Analysis", selectionName = "Feature\nSelection"), fill = switch(boxFillColouring, validation = "Validation", datasetName = "Dataset", classificationName = "Analysis", selectionName = "Feature\nSelection"))
  
  if(rowVariable != "None" || columnVariable != "None")
    performancePlot <- performancePlot + ggplot2::facet_grid(reformulate(switch(columnVariable, validation = "validation", datasetName = "dataset", classificationName = "analysis", selectionName = "selection", None = '.'), switch(rowVariable, validation = "validation", datasetName = "dataset", classificationName = "analysis", selectionName = "selection", None = '.'))) + ggplot2::theme(strip.text = ggplot2::element_text(size = fontSizes[4]))

  if(plot == TRUE)
    print(performancePlot)
  
  performancePlot
})