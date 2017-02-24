setGeneric("samplesMetricMap", function(results, ...)
{standardGeneric("samplesMetricMap")})

setMethod("samplesMetricMap", "list", 
          function(results,
                   comparison = c("classificationName", "datasetName", "selectionName", "validation"),
                   metric = c("error", "accuracy"),
                   metricColours = list(c("#0000FF", "#3F3FFF", "#7F7FFF", "#BFBFFF", "#FFFFFF"),
                                       c("#FF0000", "#FF3F3F", "#FF7F7F", "#FFBFBF", "#FFFFFF")),
                   classColours = c("blue", "red"), fontSizes = c(24, 16, 12, 12, 12),
                   mapHeight = 4, title = "Error Comparison", showLegends = TRUE, xAxisLabel = "Sample Name", showXtickLabels = TRUE,
                   showYtickLabels = TRUE, yAxisLabel = "Analysis", legendSize = grid::unit(1, "lines"), plot = TRUE)
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")  
  if(!requireNamespace("gridExtra", quietly = TRUE))
    stop("The package 'gridExtra' could not be found. Please install it.")       
  if(!requireNamespace("gtable", quietly = TRUE))
    stop("The package 'gtable' could not be found. Please install it.")
            
  comparison <- match.arg(comparison)
  metric <- match.arg(metric)
  metricText <- switch(metric, error = "Error", accuracy = "Accuracy")
            
  nColours <- if(is.list(metricColours)) length(metricColours[[1]]) else length(metricColours)
  metricBinEnds <- seq(0, 1, 1/nColours)
  knownClasses <- actualClasses(results[[1]])
  
  metricValues <- lapply(results, function(result)
  {
    resultTable <- do.call(rbind, predictions(result))
    sampleMetricValues <- by(resultTable, resultTable[, "sample"],
                          function(samplePredictions)
                          {
                            predictedClasses <- samplePredictions[, "label"]
                            actualClasses <- result@actualClasses[samplePredictions[1, "sample"]]
                            if(metric == "error")
                              sum(predictedClasses != actualClasses)
                            else
                              sum(predictedClasses == actualClasses)
                          })
                         
    cut(sampleMetricValues / table(resultTable[, "sample"]), metricBinEnds, include.lowest = TRUE)
  })
  classedMetricValues <- lapply(metricValues, function(metricSet)
  {
    metricSet <- factor(paste(knownClasses, metricSet, sep = ','),
                       levels = c(t(outer(levels(knownClasses), levels(metricSet), paste, sep = ','))))
  })
  
  ordering <- order(knownClasses)
  knownClasses <- knownClasses[ordering]
  metricValues <- lapply(metricValues, function(resultmetricValues) resultmetricValues[ordering])
  classedMetricValues <- lapply(classedMetricValues, function(resultmetricValues) resultmetricValues[ordering])
  compareFactor <- switch(comparison, classificationName = sapply(results, function(result) result@classificationName),
                                      datasetName = sapply(results, function(result) result@datasetName),
                                      selectionName = sapply(results, function(result) result@selectResult@selectionName),
                                      validation = sapply(results, function(result) .validationText(result)))
  
  plotData <- data.frame(name = factor(rep(sampleNames(results[[1]])[ordering], length(results)), levels = sampleNames(results[[1]])[ordering]),
                         type = factor(rep(compareFactor, sapply(metricValues, length)), levels = rev(compareFactor)),
                         class = rep(knownClasses, length(results)),
                         Metric = unlist(metricValues))
  
  originalLegends <- showLegends                       
  originalmetricColours <- metricColours
  showLegends <- FALSE
  if(is.list(metricColours))
  {
    metricColours <- unlist(metricColours)
    plotData[, "Metric"] <- unlist(classedMetricValues)
  }                         
  
  classData <- data.frame(Class = knownClasses)
  classesPlot <- ggplot2::ggplot(classData, ggplot2::aes(1:length(knownClasses), factor(1)), environment = environment()) +
    ggplot2::scale_fill_manual(values = classColours) + ggplot2::geom_tile(ggplot2::aes(fill = Class)) +
    ggplot2::scale_x_discrete(expand = c(0, 0), breaks = NULL, limits = c(1, length(knownClasses))) +
    ggplot2::scale_y_discrete(expand = c(0, 0), breaks = NULL) +
    ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(0.2, 0, -1, 0), "lines"),
                                                   legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                   legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                   legend.position = ifelse(showLegends, "right", "none"),
                                                   legend.key.size = legendSize)
                                                           
  metricPlot <- ggplot2::ggplot(plotData, ggplot2::aes(name, type)) + ggplot2::geom_tile(ggplot2::aes(fill = Metric)) +
    ggplot2::scale_fill_manual(values = metricColours, drop = FALSE) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme_bw() +
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                   axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(angle = 45, hjust = 1, size = fontSizes[3]) else ggplot2::element_blank(),
                   axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3]) else ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(size = fontSizes[2]),
                   axis.title.y = ggplot2::element_text(size = fontSizes[2]),
                   plot.margin = grid::unit(c(0, 1, 1, 1), "lines"),
                   legend.title = ggplot2::element_text(size = fontSizes[4]),
                   legend.text = ggplot2::element_text(size = fontSizes[5]),
                   legend.position = ifelse(showLegends, "right", "none"),
                   legend.key.size = legendSize) + ggplot2::labs(x = xAxisLabel, y = yAxisLabel)
  
  classGrob <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(classesPlot))
  metricGrob <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(metricPlot))
  commonWidth <- grid::unit.pmax(classGrob[["widths"]], metricGrob[["widths"]])
  classGrob[["widths"]] <- commonWidth
  metricGrob[["widths"]] <- commonWidth
  
  if(originalLegends == TRUE)
  {
    showLegends <- TRUE
    metricColours <- originalmetricColours
    
    classesPlot <- ggplot2::ggplot(classData, ggplot2::aes(1:length(knownClasses), factor(1)), environment = environment()) +
      ggplot2::scale_fill_manual(values = classColours) + ggplot2::geom_tile(ggplot2::aes(fill = Class)) +
      ggplot2::scale_x_discrete(expand = c(0, 0), breaks = NULL, limits = c(1, length(knownClasses))) +
      ggplot2::scale_y_discrete(expand = c(0, 0), breaks = NULL) +
      ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(0.2, 0, -1, 0), "lines"),
                                                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                     legend.position = ifelse(showLegends, "right", "none"),
                                                     legend.key.size = legendSize)
    
    classGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(classesPlot))            
    plotData[, "Metric"] <- unlist(metricValues)
    metricPlot <- ggplot2::ggplot(plotData, ggplot2::aes(name, type)) + ggplot2::geom_tile(ggplot2::aes(fill = Metric)) +
      ggplot2::scale_fill_manual(name = paste(levels(knownClasses)[1], metricText),
                                 values = if(is.list(metricColours)) metricColours[[1]] else metricColours, drop = FALSE) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme_bw() +
      ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                     axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(angle = 45, hjust = 1, size = fontSizes[3]) else ggplot2::element_blank(),
                     axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3]) else ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_text(size = fontSizes[2]),
                     axis.title.y = ggplot2::element_text(size = fontSizes[2]),
                     plot.margin = grid::unit(c(0, 1, 1, 1), "lines"),
                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                     legend.position = ifelse(showLegends, "right", "none"),
                     legend.key.size = legendSize) + ggplot2::labs(x = xAxisLabel, y = yAxisLabel)
    
    metricGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(metricPlot))
    commonWidth <- grid::unit.pmax(classGrobUnused[["widths"]], metricGrobUnused[["widths"]])                   
    metricGrobUnused[["widths"]] <- commonWidth  
    classGrobUnused[["widths"]] <- commonWidth  
    classLegend <- classGrobUnused[["grobs"]][[which(sapply(classGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
    if(showLegends == TRUE)    
      firstLegend <- metricGrobUnused[["grobs"]][[which(sapply(metricGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
    
    metricPlot <- ggplot2::ggplot(plotData, ggplot2::aes(name, type)) + ggplot2::geom_tile(ggplot2::aes(fill = Metric)) +
      ggplot2::scale_fill_manual(name = paste(levels(knownClasses)[2], metricText),
                                 values = if(is.list(metricColours)) metricColours[[2]] else metricColours, drop = FALSE) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme_bw() +
      ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                     axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(angle = 45, hjust = 1, size = fontSizes[3]) else ggplot2::element_blank(),
                     axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3]) else ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_text(size = fontSizes[2]),
                     axis.title.y = ggplot2::element_text(size = fontSizes[2]),
                     plot.margin = grid::unit(c(0, 1, 1, 1), "lines"),
                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                     legend.position = ifelse(showLegends, "right", "none"),
                     legend.key.size = legendSize) + ggplot2::labs(x = xAxisLabel, y = yAxisLabel)
    
    metricGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(metricPlot)) 
    metricGrobUnused[["widths"]] <- commonWidth  
    secondLegend <- metricGrobUnused[["grobs"]][[which(sapply(metricGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
  }
  
  if(showLegends == TRUE)
  {
    legendWidth <- max(sum(classLegend[["widths"]]), sum(firstLegend[["widths"]]))
    legendHeight <- firstLegend[["heights"]][3]
    widths <- unit.c(unit(1, "npc") - legendWidth, legendWidth)
    heights <- unit.c(unit(1 / (mapHeight + 1), "npc"), legendHeight)
    if(is.list(metricColours))
    {
      heights <- unit.c(heights, legendHeight)
      heights <- unit.c(heights, unit(1, "npc") - 2 * legendHeight - unit(1 / (mapHeight + 1), "npc"))
    } else
    {
      heights <- unit.c(heights, unit(1, "npc") - legendHeight - unit(1 / (mapHeight + 1), "npc"))
    }
  }
  else
  {
    widths <- unit(1, "npc")
    heights <- unit.c(unit(1 / (mapHeight + 1), "npc"), unit(1, "npc") - unit(1 / (mapHeight + 1), "npc"))
  }
  
  grobTable <- gtable::gtable(widths, heights)
  grobTable <- gtable::gtable_add_grob(grobTable, classGrob, 1, 1)
  if(showLegends == TRUE)
  {
    if(is.list(metricColours))
      grobTable <- gtable::gtable_add_grob(grobTable, metricGrob, 2, 1, 4, 1)
    else
      grobTable <- gtable::gtable_add_grob(grobTable, metricGrob, 2, 1, 3, 1)
  } else
  {
    grobTable <- gtable::gtable_add_grob(grobTable, metricGrob, 2, 1)
  }
  if(showLegends == TRUE)
  {  
    grobTable <- gtable::gtable_add_grob(grobTable, classLegend, 1, 2)
    grobTable <- gtable::gtable_add_grob(grobTable, firstLegend, 2, 2)
  }
  if(showLegends == TRUE && is.list(metricColours))
  {
    grobTable <- gtable::gtable_add_grob(grobTable, secondLegend, 3, 2)
    grobTable <- gtable::gtable_add_grob(grobTable, grob(), 4, 2)
  }
  wholePlot <- gridExtra::arrangeGrob(grobTable, top = grid::textGrob(title, gp = grid::gpar(fontsize = fontSizes[1])))
  
  if(plot == TRUE)               
    grid::grid.draw(wholePlot)
  wholePlot
})