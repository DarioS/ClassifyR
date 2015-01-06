setGeneric("performancePlot", function(results, ...)
{standardGeneric("performancePlot")})

setMethod("performancePlot", "list", 
          function(results,
                   aggregate = FALSE,
                   xVariable = c("classificationName", "datasetName", "validation"),
                   performanceName = NULL,
                   boxFillColouring = c("classificationName", "datasetName", "validation", "None"),
                   boxFillColours = NULL,
                   boxLineColouring = c("classificationName", "datasetName", "validation", "None"),
                   boxLineColours = NULL,
                   yMax = 1, fontSizes = c(24, 16, 12), title = NULL,
                   xLabel = "Classification", yLabel = performanceName,
                   margin = grid::unit(c(0, 1, 1, 0), "lines"), rotate90 = FALSE, plot = TRUE)
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")             
  if(!requireNamespace("scales", quietly = TRUE))
    stop("The package 'scales' could not be found. Please install it.")
  
  xVariable <- match.arg(xVariable)
  boxFillColouring <- match.arg(boxFillColouring)
  boxLineColouring <- match.arg(boxLineColouring)
  if(is.null(performanceName)) stop("Please specify a performance measure to plot.")
  
  performances <- lapply(results, function(result)
                         {
                           if(!performanceName %in% names(result@performance))
                             stop("Performance measure not calculated for ", result@classificationName)
                           if(aggregate == TRUE) mean(result@performance[[performanceName]])
                           else result@performance[[performanceName]]
                         })
  performanceLengths <- sapply(performances, length)
  
  plotData <- data.frame(dataset = rep(sapply(results, function(result) result@datasetName), performanceLengths),
                         analysis = rep(sapply(results, function(result) result@classificationName), performanceLengths),
                         validation = rep(sapply(results, function(result) .validationText(result)), performanceLengths),
                         performance = unlist(performances))  
  
  if(boxFillColouring != "None")
    if(is.null(boxFillColours)) boxFillColours <- scales::hue_pal()(switch(boxFillColouring, validation = length(unique(plotData[, "validation"])), datasetName = length(unique(plotData[, "dataset"])), classificationName = length(unique(plotData[, "analysis"]))))
  if(boxLineColouring != "None")
    if(is.null(boxLineColours)) boxLineColours <- scales::hue_pal(direction = -1)(switch(boxLineColouring, validation = length(unique(plotData[, "validation"])), datasetName = length(unique(plotData[, "dataset"])), classificationName = length(unique(plotData[, "analysis"]))))
  
  plotData[, "dataset"] <- factor(plotData[, "dataset"], levels = unique(sapply(results, function(result) result@datasetName)))
  plotData[, "analysis"] <- factor(plotData[, "analysis"], levels = unique(sapply(results, function(result) result@classificationName)))
  plotData[, "validation"] <- factor(plotData[, "validation"], levels = unique(plotData[, "validation"])) # Order in which validations were provided.
  
  if(rotate90 == TRUE)
  {
    switch(xVariable, validation = plotData[, "validation"] <- factor(plotData[, "validation"], levels = rev(levels(plotData[, "validation"]))),
           datasetName = plotData[, "dataset"] <- factor(plotData[, "dataset"], levels = rev(levels(plotData[, "dataset"]))),
           classificationName = plotData[, "analysis"] <- factor(plotData[, "analysis"], levels = rev(levels(plotData[, "analysis"]))))
  }

  performancePlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = switch(xVariable, validation = validation, datasetName = dataset, classificationName = analysis), y = performance,
                          fill = switch(boxFillColouring, validation = validation, datasetName = dataset, classificationName = analysis, None = NULL), colour = switch(boxLineColouring, validation = validation, datasetName = dataset, classificationName = analysis, None = NULL)), environment = environment()) +
                          ggplot2::scale_y_continuous(limits = c(0, yMax)) + ggplot2::scale_fill_manual(values = boxFillColours) + ggplot2::scale_colour_manual(values = boxLineColours) + ggplot2::xlab(xLabel) + ggplot2::ylab(yLabel) +
                          ggplot2::ggtitle(title) + ggplot2::theme(legend.position = "none", axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), plot.title = ggplot2::element_text(size = fontSizes[1]), plot.margin = margin)
  geom <- switch(xVariable,
         validation = if(length(plotData[, "validation"]) == length(unique(plotData[, "validation"]))) ggplot2::geom_bar(stat = "identity") else ggplot2::geom_boxplot(),
         datasetName = if(length(plotData[, "dataset"]) == length(unique(plotData[, "dataset"]))) ggplot2::geom_bar(stat = "identity") else ggplot2::geom_boxplot(),
         classificationName = if(length(plotData[, "analysis"]) == length(unique(plotData[, "analysis"]))) ggplot2::geom_bar(stat = "identity") else ggplot2::geom_boxplot())
  performancePlot <- performancePlot + geom
  if(rotate90 == TRUE)
    performancePlot <- performancePlot + ggplot2::coord_flip()
  
  if(plot == TRUE)
    print(performancePlot)
  
  performancePlot
})