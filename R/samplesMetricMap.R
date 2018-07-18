setGeneric("samplesMetricMap", function(results, ...)
{standardGeneric("samplesMetricMap")})

setMethod("samplesMetricMap", "list", 
          function(results,
                   comparison = c("classificationName", "datasetName", "selectionName", "validation"),
                   metric = c("error", "accuracy"),
                   featureValues = NULL, featureName = NULL,
                   metricColours = list(c("#3F48CC", "#6F75D8", "#9FA3E5", "#CFD1F2", "#FFFFFF"),
                                        c("#880015", "#A53F4F", "#C37F8A", "#E1BFC4", "#FFFFFF")),
                   classColours = c("#3F48CC", "#880015"), groupColours = c("darkgreen", "yellow2"),
                   fontSizes = c(24, 16, 12, 12, 12),
                   mapHeight = 4, title = "Error Comparison", showLegends = TRUE, xAxisLabel = "Sample Name", showXtickLabels = TRUE,
                   yAxisLabel = "Analysis", showYtickLabels = TRUE, legendSize = grid::unit(1, "lines"), plot = TRUE)
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
  metricID <- switch(metric, error = "Sample-wise Error Rate", accuracy = "Sample-wise Accuracy")
  allCalculated <- all(sapply(results, function(result) metricID %in% names(performance(result))))
  if(!allCalculated)
    stop("One or more classification results lack the calculated sample-specific metric.")
  if(!is.null(featureValues) && is.null(featureName))
    stop("featureValues is specified by featureNames isn't. Specify both.")
  
  nColours <- if(is.list(metricColours)) length(metricColours[[1]]) else length(metricColours)
  metricBinEnds <- seq(0, 1, 1/nColours)
  knownClasses <- actualClasses(results[[1]])

  metricValues <- lapply(results, function(result)
  {
    sampleMetricValues <- result@performance[[metricID]]
    cut(sampleMetricValues, metricBinEnds, include.lowest = TRUE)
  })
  classedMetricValues <- lapply(metricValues, function(metricSet)
  {
    metricSet <- factor(paste(knownClasses, metricSet, sep = ','),
                       levels = c(t(outer(levels(knownClasses), levels(metricSet), paste, sep = ','))))
  })
  
  if(is.null(featureValues))
  {    
    ordering <- order(knownClasses)
  } else {
    featureValues <- featureValues[match(sampleNames(results[[1]]), names(featureValues))]
    ordering <- order(knownClasses, featureValues)
  }
  
  knownClasses <- knownClasses[ordering]
  if(!is.null(featureValues))
    featureValues <- featureValues[ordering]
  
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
    ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(0.2, 0, 0.01, 0), "npc"),
                                                   legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                   legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                   legend.position = ifelse(showLegends, "right", "none"),
                                                   legend.key.size = legendSize)
  
  if(!is.null(featureValues))
  {
    if(is.factor(featureValues))
    {
      featureValuesData <- data.frame(Group = featureValues)
      featureValuesPlot <- ggplot2::ggplot(featureValuesData, ggplot2::aes(1:length(featureValues), factor(1)), environment = environment()) +
      ggplot2::scale_fill_manual(name = featureName, values = groupColours) + ggplot2::geom_tile(ggplot2::aes(fill = Group)) +
      ggplot2::scale_x_discrete(expand = c(0, 0), breaks = NULL, limits = c(1, length(featureValues))) +
      ggplot2::scale_y_discrete(expand = c(0, 0), breaks = NULL) +
      ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "npc"),
                                                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                     legend.position = ifelse(showLegends, "right", "none"),
                                                     legend.key.size = legendSize)
    } else # Numeric data about the samplles.
    {
      featureValuesData <- data.frame(Class = knownClasses, measurements = featureValues)
      featureValuesPlot <- ggplot2::ggplot(featureValuesData, environment = environment()) +
      ggplot2::geom_point(ggplot2::aes(x = 1:length(featureValues), y = measurements, colour = Class)) +
      ggplot2::scale_colour_manual(values = classColours) +
      ggplot2::scale_x_discrete(breaks = NULL, limits = c(1, length(featureValues))) +
      ggplot2::scale_y_continuous(breaks = c(min(featureValues), max(featureValues))) +
      ggplot2::labs(x = featureName, y = '') + ggplot2::theme(plot.margin = grid::unit(c(0, 0, 1, 0), "lines"),
                                                     axis.text = ggplot2::element_text(colour = "black"),
                                                     axis.title.x = ggplot2::element_text(size = fontSizes[4]),          
                                                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                     legend.position = ifelse(showLegends, "right", "none"),
                                                     legend.key.size = legendSize, panel.background = ggplot2::element_blank(),
                                                     panel.border = ggplot2::element_rect(colour = "black", fill = "transparent")
                                                     )      
    }
  }
                                                   
  metricPlot <- ggplot2::ggplot(plotData, ggplot2::aes(name, type)) + ggplot2::geom_tile(ggplot2::aes(fill = Metric)) +
    ggplot2::scale_fill_manual(values = metricColours, na.value = "grey", drop = FALSE) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme_bw() +
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                   axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(angle = 45, hjust = 1, size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                   axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(size = fontSizes[2]),
                   axis.title.y = ggplot2::element_text(size = fontSizes[2]),
                   plot.margin = grid::unit(c(0, 1, 1, 1), "lines"),
                   legend.title = ggplot2::element_text(size = fontSizes[4]),
                   legend.text = ggplot2::element_text(size = fontSizes[5]),
                   legend.position = ifelse(showLegends, "right", "none"),
                   legend.key.size = legendSize) + ggplot2::labs(x = xAxisLabel, y = yAxisLabel)

  classGrob <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(classesPlot))
  metricGrob <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(metricPlot))
  if(!is.null(featureValues))
    featureValuesGrob <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(featureValuesPlot))
    
  if(!is.null(featureValues))
    commonWidth <- grid::unit.pmax(classGrob[["widths"]], metricGrob[["widths"]], featureValuesGrob[["widths"]])
  else
    commonWidth <- grid::unit.pmax(classGrob[["widths"]], metricGrob[["widths"]])
  classGrob[["widths"]] <- commonWidth
  metricGrob[["widths"]] <- commonWidth
  if(!is.null(featureValues))
    featureValuesGrob[["widths"]] <- commonWidth
  
  if(originalLegends == TRUE)
  {
    showLegends <- TRUE
    metricColours <- originalmetricColours

    classesPlot <- ggplot2::ggplot(classData, ggplot2::aes(1:length(knownClasses), factor(1)), environment = environment()) +
      ggplot2::scale_fill_manual(values = classColours) + ggplot2::geom_tile(ggplot2::aes(fill = Class)) +
      ggplot2::scale_x_discrete(expand = c(0, 0), breaks = NULL, limits = c(1, length(knownClasses))) +
      ggplot2::scale_y_discrete(expand = c(0, 0), breaks = NULL) +
      ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(0.2, 0, 0.01, 0), "lines"),
                                                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                     legend.position = ifelse(showLegends, "right", "none"),
                                                     legend.key.size = legendSize)
    classGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(classesPlot)) 
    
    if(!is.null(featureValues) && is.factor(featureValues))
    {
      featureValuesPlot <- ggplot2::ggplot(featureValuesData, ggplot2::aes(1:length(featureValues), factor(1)), environment = environment()) +
        ggplot2::scale_fill_manual(name = featureName, values = groupColours) + ggplot2::geom_tile(ggplot2::aes(fill = Group)) +
        ggplot2::scale_x_discrete(expand = c(0, 0), breaks = NULL, limits = c(1, length(featureValues))) +
        ggplot2::scale_y_discrete(expand = c(0, 0), breaks = NULL) +
        ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "lines"),
                                                       legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                       legend.text = ggplot2::element_text(size = fontSizes[5]),
                                                       legend.position = ifelse(showLegends, "right", "none"),
                                                       legend.key.size = legendSize)
      featureValuesGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(featureValuesPlot))     
    } else {featureValuesGrobUnused <- grob()}
    
    plotData[, "Metric"] <- unlist(metricValues)
    classLegend <- NULL
    if(is.list(metricColours))
      classLegend <- paste(levels(knownClasses)[1], NULL)
    metricPlot <- ggplot2::ggplot(plotData, ggplot2::aes(name, type)) + ggplot2::geom_tile(ggplot2::aes(fill = Metric)) +
      ggplot2::scale_fill_manual(name = paste(classLegend, metricText, sep = ''),
                                 values = if(is.list(metricColours)) metricColours[[1]] else metricColours, na.value = "grey", drop = FALSE) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme_bw() +
      ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                     axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(angle = 45, hjust = 1, size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                     axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_text(size = fontSizes[2]),
                     axis.title.y = ggplot2::element_text(size = fontSizes[2]),
                     plot.margin = grid::unit(c(0, 1, 0, 1), "lines"),
                     legend.title = ggplot2::element_text(size = fontSizes[4]),
                     legend.text = ggplot2::element_text(size = fontSizes[5]),
                     legend.position = ifelse(showLegends, "right", "none"),
                     legend.key.size = legendSize) + ggplot2::labs(x = xAxisLabel, y = yAxisLabel)
    
    metricGrobUnused <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(metricPlot))
    if(!is.null(featureValues) && is.factor(featureValues))
      commonWidth <- grid::unit.pmax(classGrobUnused[["widths"]], metricGrobUnused[["widths"]], featureValuesGrobUnused[["widths"]])                   
    else
      commonWidth <- grid::unit.pmax(classGrobUnused[["widths"]], metricGrobUnused[["widths"]])                   
    metricGrobUnused[["widths"]] <- commonWidth  
    classGrobUnused[["widths"]] <- commonWidth
    if(!is.null(featureValues) && is.factor(featureValues))
      featureValuesGrobUnused[["widths"]] <- commonWidth
    classLegend <- classGrobUnused[["grobs"]][[which(sapply(classGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
    if(!is.null(featureValues) && is.factor(featureValues))
      featureValuesLegend <- featureValuesGrobUnused[["grobs"]][[which(sapply(featureValuesGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
    else
      featureValuesLegend <- grob()
    if(showLegends == TRUE)    
      firstLegend <- metricGrobUnused[["grobs"]][[which(sapply(metricGrobUnused[["grobs"]], function(grob) grob[["name"]]) == "guide-box")]]
    
    metricPlot <- ggplot2::ggplot(plotData, ggplot2::aes(name, type)) + ggplot2::geom_tile(ggplot2::aes(fill = Metric)) +
      ggplot2::scale_fill_manual(name = paste(levels(knownClasses)[2], metricText),
                                 values = if(is.list(metricColours)) metricColours[[2]] else metricColours, na.value = "grey", drop = FALSE) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme_bw() +
      ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                     axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(angle = 45, hjust = 1, size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                     axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
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
    classesHeight <- grid::unit(1 / (mapHeight + 1), "npc")
    if(!is.null(featureValues) && is.factor(featureValues))
      legendWidth <- max(sum(classLegend[["widths"]]), sum(firstLegend[["widths"]]), sum(featureValuesLegend[["widths"]]))
    else
      legendWidth <- max(sum(classLegend[["widths"]]), sum(firstLegend[["widths"]]))
    featureValuesHeight <- grid::unit(0, "cm")
    if(!is.null(featureValues))
    {
      if(is.factor(featureValues))
        featureValuesHeight <- featureValuesLegend[["heights"]][3]
      else
        featureValuesHeight <- classesHeight
    }
    
    if(is.list(metricColours))
      legendHeight <- (grid::unit(1, "npc") - featureValuesHeight) * (mapHeight / (mapHeight + 1) / 2)
    else
      legendHeight <- (grid::unit(1, "npc") - featureValuesHeight) * (mapHeight / (mapHeight + 1))
    
    widths <- grid::unit.c(unit(1, "npc") - legendWidth, legendWidth)
    heights <- grid::unit.c(unit(1 / (mapHeight + 1), "npc"), featureValuesHeight, legendHeight)
    if(is.list(metricColours))
    {
      heights <- grid::unit.c(heights, legendHeight)
    } else # Greyscale legend.
    {
      heights <- grid::unit.c(heights, unit(1, "npc") - legendHeight - featureValuesHeight - grid::unit(1 / (mapHeight + 1), "npc"))
    }
  }
  else
  {
    widths <- grid::unit(1, "npc")
    heights <- grid::unit.c(grid::unit(1 / (mapHeight + 1), "npc"), featureValuesHeight, grid::unit(1, "npc") - featureValuesHeight - grid::unit(1 / (mapHeight + 1), "npc"))
  }
  
  classLegend[["vp"]][["valid.just"]] <- c(0.7, 0.5)
  if(!is.null(featureValues) && is.factor(featureValues))
    featureValuesLegend[["vp"]][["valid.just"]] <- c(0.7, 0.33)

  grobTable <- gtable::gtable(widths, heights)
  grobTable <- gtable::gtable_add_grob(grobTable, classGrob, 1, 1)
  if(!is.null(featureValues))
    grobTable <- gtable::gtable_add_grob(grobTable, featureValuesGrob, 2, 1)
  else
    grobTable <- gtable::gtable_add_grob(grobTable, grob(), 2, 1)
  if(showLegends == TRUE)
  {
    if(is.list(metricColours))
      grobTable <- gtable::gtable_add_grob(grobTable, metricGrob, 3, 1, 4, 1)
    else
      grobTable <- gtable::gtable_add_grob(grobTable, metricGrob, 3, 1, 3, 1)
  } else
  {
    grobTable <- gtable::gtable_add_grob(grobTable, metricGrob, 3, 1)
  }
  if(showLegends == TRUE)
  {  
    grobTable <- gtable::gtable_add_grob(grobTable, classLegend, 1, 2)
    grobTable <- gtable::gtable_add_grob(grobTable, featureValuesLegend, 2, 2)
  }
  if(showLegends == TRUE && !is.list(metricColours))
  {
    firstLegend[["vp"]][["valid.just"]] <- c(0.62, 0.5)
    grobTable <- gtable::gtable_add_grob(grobTable, firstLegend, 3, 2)
  }
  if(showLegends == TRUE && is.list(metricColours))
  {
    firstLegend[["vp"]][["valid.just"]] <- c(0.62, 0.5)
    secondLegend[["vp"]][["valid.just"]] <- c(0.62, 0.4)
    grobTable <- gtable::gtable_add_grob(grobTable, firstLegend, 3, 2)
    grobTable <- gtable::gtable_add_grob(grobTable, secondLegend, 4, 2)
  }
  wholePlot <- gridExtra::arrangeGrob(grobTable, top = grid::textGrob(title, vjust = 0.5, gp = grid::gpar(fontsize = fontSizes[1])))

  if(plot == TRUE)               
    grid::grid.draw(wholePlot)
  wholePlot
})