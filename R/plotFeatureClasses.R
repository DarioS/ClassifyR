setGeneric("plotFeatureClasses", function(expression, ...)
           {standardGeneric("plotFeatureClasses")})

setMethod("plotFeatureClasses", "matrix", 
          function(expression, classes, ...)
{
  groupsTable <- data.frame(class = classes)
  features <- rownames(expression)
  rownames(expression) <- NULL
  rownames(groupsTable) <- colnames(expression)
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features
  plotFeatureClasses(exprSet, ...)
})

setMethod("plotFeatureClasses", "ExpressionSet", 
          function(expression, rows, whichPlots = c("both", "density", "stripchart"),
                   xAxisLabel = expression(log[2](expression)), expressionLimits = c(2, 16),
                   yAxisLabels = c("Density", "Classes"),
                   showXtickLabels = TRUE, showYtickLabels = TRUE, xLabelPositions = "auto",
                   yLabelPositions = "auto", fontSizes = c(24, 16, 12, 12, 12),
                   colours = c("blue", "red"), plot = TRUE)
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")  
  if(!requireNamespace("grid", quietly = TRUE))
    stop("The package 'grid' could not be found. Please install it.")
  if(!requireNamespace("gridExtra", quietly = TRUE))
    stop("The package 'gridExtra' could not be found. Please install it.")     

  whichPlots <- match.arg(whichPlots)
  classes <- pData(expression)[, "class"]
  features <- featureNames(expression)
  expression <- exprs(expression)
  if(xLabelPositions[1] == "auto")
    xLabelPositions <- ggplot2::waiver()
  if(yLabelPositions[1] == "auto")
    yLabelPositions <- ggplot2::waiver()
  
  invisible(lapply(rows, function(featureRow)
  {
    plotData <- data.frame(expr = expression[featureRow, ], classes)
    if(whichPlots %in% c("both", "density"))
    {
      densPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = expr, colour = classes)) +
        ggplot2::stat_density(ggplot2::aes(y = ..density..), geom = "path", position = "identity", size = 1) +
        ggplot2::scale_colour_manual("Class", values = c("red", "blue")) + ggplot2::coord_cartesian(xlim = expressionLimits) +
        ggplot2::scale_x_continuous(breaks = xLabelPositions) + ggplot2::scale_y_continuous(breaks = yLabelPositions)
    }

    if(whichPlots %in% c("both", "stripchart"))
    {
      stripPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = classes, y = expr)) +
                   ggplot2::geom_dotplot(dotsize = 0.5, binaxis = "y", stackdir = "center", position = "dodge", ggplot2::aes(colour = classes)) +
                   ggplot2::scale_colour_manual("Class", values = c("red", "blue")) + ggplot2::xlab(yAxisLabels[2]) + ggplot2::ylab(xAxisLabel) + ggplot2::scale_y_continuous(limits = expressionLimits) +
                   ggplot2::theme(plot.title = ggplot2::element_text(size = fontSizes[1]),
                                  axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                                  axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                                  axis.ticks.x = if(showXtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                                  axis.title = ggplot2::element_text(size = fontSizes[2], colour = "black"), legend.title = ggplot2::element_text(size = fontSizes[4]),
                                  legend.text = ggplot2::element_text(size = fontSizes[5]))
      stripPlot <- stripPlot + ggplot2::coord_flip()
    }
    if(is.numeric(rows)) featureName <- features[featureRow] else featureName <-featureRow
    if(whichPlots == "both")
    {
      densPlot <- densPlot + ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
                                            axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                                            axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                                            axis.ticks.x = if(showXtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                                            axis.ticks.y = if(showYtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                                            legend.title = ggplot2::element_text(size = fontSizes[4]),
                                            legend.text = ggplot2::element_text(size = fontSizes[5])) +
                             ggplot2::labs(x = NULL, y = yAxisLabels[1])
      densityGraphicsTable <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(densPlot))
      stripGraphicsTable <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(stripPlot))
      maxWidth = grid::unit.pmax(densityGraphicsTable[["widths"]][2:3], stripGraphicsTable[["widths"]][2:3])
      densityGraphicsTable[["widths"]][2:3] <- maxWidth
      stripGraphicsTable[["widths"]][2:3] <- maxWidth
      
      bothGraphics <- gridExtra::arrangeGrob(densityGraphicsTable, stripGraphicsTable, nrow = 2,
                                              top = grid::textGrob(featureName, gp = grid::gpar(fontsize = fontSizes[1]), vjust = 1))
      if(plot == TRUE)
      {
        grid::grid.draw(bothGraphics)
        if(featureRow != rows[length(rows)])
          grid::grid.newpage()
      }
      bothGraphics
    } else if(whichPlots == "density")
    {
      densPlot <- densPlot + ggplot2::xlab(xAxisLabel) +  ggplot2::ylab(yAxisLabels[1]) +
                  ggplot2::theme(plot.title = ggplot2::element_text(size = fontSizes[1]),
                                 axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                                 axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                                 axis.ticks.x = if(showXtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                                 axis.ticks.y = if(showYtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                                 axis.title = ggplot2::element_text(size = fontSizes[2], colour = "black"), legend.title = ggplot2::element_text(size = fontSizes[4]),
                                 legend.text = ggplot2::element_text(size = fontSizes[5])) + ggplot2::ggtitle(featureName)
                  
      if(plot == TRUE)
        print(densPlot)
      densPlot
    } else {
      stripPlot <- stripPlot + ggplot2::ggtitle(featureName)
      if(plot == TRUE)
        print(stripPlot)
      stripPlot
    }
  }))
})
