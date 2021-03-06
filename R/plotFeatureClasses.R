setGeneric("plotFeatureClasses", function(measurements, ...)
           {standardGeneric("plotFeatureClasses")})

setMethod("plotFeatureClasses", "matrix", function(measurements, classes, targets, ...)
{
  if(missing(targets))
    stop("'targets' must be specified.")
  
  if(!"Pairs" %in% class(targets))
    flippedMatrix <- t(measurements[targets, , drop = FALSE])
  else
    flippedMatrix <- t(measurements[union(S4Vectors::first(targets), S4Vectors::second(targets)), , drop = FALSE])
  plotFeatureClasses(DataFrame(flippedMatrix, check.names = FALSE), classes, targets, ...)
})

setMethod("plotFeatureClasses", "DataFrame", function(measurements, classes, targets, groupBy = NULL,
                                groupingName = NULL, whichNumericFeaturePlots = c("both", "density", "stripchart"),
                                measurementLimits = NULL, lineWidth = 1, dotBinWidth = 1,
                                xAxisLabel = NULL, yAxisLabels = c("Density", "Classes"),
                                showXtickLabels = TRUE, showYtickLabels = TRUE,
                                xLabelPositions = "auto", yLabelPositions = "auto",
                                fontSizes = c(24, 16, 12, 12, 12),
                                colours = c("#3F48CC", "#880015"),
                                showDatasetName = TRUE, plot = TRUE)
{
  if(missing(targets))
    stop("'targets' must be specified.")

  if(is.character(groupBy))
  {
    groupingName <- groupBy
    groupBy <- measurements[, groupBy]
    levelsOrder <- levels(groupBy)
    groupBy <- list(legends = factor(groupBy, levels = levelsOrder),
                    facets = factor(paste(groupingName, "is", groupBy), levels = paste(groupingName, "is", levelsOrder)))
  }

  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  classes <- splitDataset[["classes"]]

  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")
  if(!requireNamespace("cowplot", quietly = TRUE))
    stop("The package 'cowplot' could not be found. Please install it.")  
  if(!requireNamespace("grid", quietly = TRUE))
    stop("The package 'grid' could not be found. Please install it.")
  if(!requireNamespace("gridExtra", quietly = TRUE))
    stop("The package 'gridExtra' could not be found. Please install it.")

  ggplot2::theme_set(ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA), panel.grid.major = ggplot2::element_line(colour = "grey", linetype = "dashed")))
  whichNumericFeaturePlots <- match.arg(whichNumericFeaturePlots)
  if(xLabelPositions[1] == "auto")
    xLabelPositions <- ggplot2::waiver()
  if(yLabelPositions[1] == "auto")
    yLabelPositions <- ggplot2::waiver()

  if(!is(targets, "DataFrame"))
  {
    if(!"Pairs" %in% class(targets))
      measurements <- measurements[targets]
    else
      measurements <- measurements[union(S4Vectors::first(targets), S4Vectors::second(targets))]
  }

  if(!"Pairs" %in% class(targets))
  {
    invisible(lapply(seq_along(measurements), function(columnIndex) # Plot a single feature at a time.
    {
      plotData <- data.frame(measurement = measurements[, columnIndex], class = classes)
      if(!is.null(groupBy))
      {
        plotData[, "legends grouping"] <- groupBy[["legends"]]
        plotData[, "facets grouping"] <- groupBy[["facets"]]
      }

      if(is.null(S4Vectors::mcols(measurements)))
      {
        featureText <- colnames(measurements)[columnIndex]
      } else {
        featureText <- S4Vectors::mcols(measurements)[columnIndex, "feature"]
        if(showDatasetName == TRUE)
        {
          featureText <- paste(featureText, paste('(', S4Vectors::mcols(measurements)[columnIndex, "dataset"], ')', sep = ''))
        }
      }
    
      if(is.numeric(plotData[, "measurement"]))
      {
        if(whichNumericFeaturePlots %in% c("both", "density"))
        {
          densPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = measurement, colour = class)) +
            ggplot2::stat_density(ggplot2::aes(y = ..density..), geom = "path", position = "identity", size = lineWidth) +
            ggplot2::scale_colour_manual("Class", values = colours) + ggplot2::coord_cartesian(xlim = measurementLimits) +
            ggplot2::scale_x_continuous(breaks = xLabelPositions) + ggplot2::scale_y_continuous(breaks = yLabelPositions)

          if(!is.null(groupBy))
          {
            densPlot <- densPlot + ggplot2::aes(linetype = `legends grouping`) + ggplot2::scale_linetype_discrete(name = groupingName)
          }
        }
      
        if(whichNumericFeaturePlots %in% c("both", "stripchart"))
        {
          yLabel <- ifelse(whichNumericFeaturePlots == "both", yAxisLabels[2], yAxisLabels[1])
          stripPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = class, y = measurement)) +
            ggplot2::geom_dotplot(binaxis = 'y', stackdir = "center", position = "dodge", ggplot2::aes(colour = class), binwidth = dotBinWidth) +
            ggplot2::scale_colour_manual("Class", values = colours) + ggplot2::xlab(yLabel) + ggplot2::ylab(xAxisLabel) + ggplot2::scale_y_continuous(limits = measurementLimits) +
            ggplot2::theme(plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5),
                           axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                           axis.ticks.x = if(showXtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                           axis.title = ggplot2::element_text(size = fontSizes[2], colour = "black"), legend.title = ggplot2::element_text(size = fontSizes[4]),
                           legend.text = ggplot2::element_text(size = fontSizes[5]))
          stripPlot <- stripPlot + ggplot2::coord_flip()
        
          if(!is.null(groupBy))
          {
            stripPlot <- stripPlot + ggplot2::facet_wrap(~ `facets grouping`, ncol = 1, strip.position = "left")
          }
        }
        
        if(whichNumericFeaturePlots == "both")
        {
          densPlot <- densPlot + ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
                                                axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                                                axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                                                axis.ticks.x = if(showXtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                                                axis.ticks.y = if(showYtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                                                legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                legend.text = ggplot2::element_text(size = fontSizes[5])) +
                                                ggplot2::labs(x = NULL, y = yAxisLabels[1])
          
          alignedPlots <- cowplot::align_plots(densPlot, stripPlot, align = 'v', axis = "lr")
          bothGraphics <- gridExtra::arrangeGrob(alignedPlots[[1]], alignedPlots[[2]], nrow = 2,
                                                 top = grid::textGrob(featureText, gp = grid::gpar(fontsize = fontSizes[1]), vjust = 1))
          if(plot == TRUE)
          {
            grid::grid.draw(bothGraphics)
            if(columnIndex != ncol(measurements))
              grid::grid.newpage()
          }
          bothGraphics
        } else if(whichNumericFeaturePlots == "density")
        {
          densPlot <- densPlot + ggplot2::xlab(xAxisLabel) +  ggplot2::ylab(yAxisLabels[1]) +
            ggplot2::theme(plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5),
                           axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                           axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                           axis.ticks.x = if(showXtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                           axis.ticks.y = if(showYtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                           axis.title = ggplot2::element_text(size = fontSizes[2], colour = "black"), legend.title = ggplot2::element_text(size = fontSizes[4]),
                           legend.text = ggplot2::element_text(size = fontSizes[5])) + ggplot2::ggtitle(featureText)
          
          if(plot == TRUE)
            print(densPlot)
          densPlot
        } else {
          stripPlot <- stripPlot + ggplot2::ggtitle(featureText)
          if(plot == TRUE)
            print(stripPlot)
          stripPlot
        }
      } else { # Plotting variable is a factor.
        barPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = measurement, fill = class)) +
                    ggplot2::geom_bar() + ggplot2::xlab(xAxisLabel) + ggplot2::ylab("Number of Samples") +
                    ggplot2::scale_fill_manual("Class", values = colours) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5),
                         axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                         axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                         axis.ticks.x = if(showXtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                         axis.ticks.y = if(showYtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                         axis.title = ggplot2::element_text(size = fontSizes[2], colour = "black"), legend.title = ggplot2::element_text(size = fontSizes[4]),
                         legend.text = ggplot2::element_text(size = fontSizes[5])) + ggplot2::ggtitle(featureText)
        if(!is.null(groupBy))
        {
          barPlot <- barPlot + ggplot2::facet_wrap(~ `facets grouping`, ncol = 1, strip.position = "left")
        }
        if(plot == TRUE)
          print(barPlot)
        barPlot
      }
    }))
  } else { # Plot parallel plots.
    invisible(lapply(targets, function(featurePair) # Plot a single feature at a time.
    {
      plotData <- data.frame(first = measurements[, S4Vectors::first(featurePair)],
                             second = measurements[, S4Vectors::second(featurePair)], class = classes)
      if(!is.null(groupBy))
        plotData[, "facets grouping"] <- groupBy[["facets"]]
      
      groupedPlotData <- tidyr::gather(dplyr::mutate(plotData, ID = 1:dplyr::n()), key, value, 1:2)
          
      pairsPlot <- ggplot2::ggplot(groupedPlotData, ggplot2::aes(key, value, group = ID, colour = class)) +
                   ggplot2::geom_line(size = lineWidth) +
                   ggplot2::scale_colour_manual("Class", values = colours) + ggplot2::coord_cartesian(ylim = measurementLimits) +
                   ggplot2::scale_x_discrete(expand = c(0.05, 0.05), labels = c(S4Vectors::first(featurePair), S4Vectors::second(featurePair))) +
                   ggplot2::scale_y_continuous(breaks = yLabelPositions) + ggplot2::ggtitle("Pairs Plot") +
                   ggplot2::theme(plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5),
                           axis.text.x = if(showXtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                           axis.text.y = if(showYtickLabels == TRUE) ggplot2::element_text(size = fontSizes[3], colour = "black") else ggplot2::element_blank(),
                           axis.ticks.x = if(showXtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                           axis.ticks.y = if(showYtickLabels == TRUE) ggplot2::element_line() else ggplot2::element_blank(),
                           axis.title = ggplot2::element_text(size = fontSizes[2], colour = "black"), legend.title = ggplot2::element_text(size = fontSizes[4]),
                           legend.text = ggplot2::element_text(size = fontSizes[5])) + ggplot2::labs(x = xAxisLabel, y = yAxisLabels[1])

      if(!is.null(groupBy))
      {
        pairsPlot <- pairsPlot + ggplot2::facet_wrap(~ `facets grouping`, ncol = 1, strip.position = "left")
      }
      if(plot == TRUE)
        print(pairsPlot)
      pairsPlot
    }))
  }
})

setMethod("plotFeatureClasses", "MultiAssayExperiment",
                                function(measurements, targets, groupBy = NULL, groupingName = NULL, showDatasetName = TRUE, ...)
{
  if(missing(targets))
    stop("'targets' must be specified by the user.")
  if(!all(targets[, 1] %in% c(names(measurements), "clinical")))
    stop("Some table names in 'targets' are not assay names in 'measurements' or \"clinical\".")  
                                
  assaysTargets <- targets[targets[, 1] != "clinical", ]
  sampleInfoTargets <- targets[targets[, 1] == "clinical", ]
  measurements <- measurements[assaysTargets[, 2], , assaysTargets[, 1]]
  classes <- MultiAssayExperiment::colData(measurements)[, "class"]

  if(!is.null(groupBy))
  {
    if(is.null(groupingName))
      groupingName <- groupBy[2]
    groupingTable <- groupBy[1]
    if(groupingTable == "clinical")
    {
      groupBy <- MultiAssayExperiment::colData(measurements)[, groupBy[2]]
    } else { # One of the omics tables.
      groupBy <- measurements[groupBy[2], , groupingTable]
      if(showDatasetName == TRUE)
        groupingName <- paste(groupingName, groupingTable)
    }
    levelsOrder <- levels(groupBy)
    groupBy <- list(legends = factor(groupBy, levels = levelsOrder),
                    facets = {groupText <- paste(groupingName, "is", groupBy)
                             factor(groupText, levels = paste(groupingName, "is", levelsOrder))}
                    )
  }

  MultiAssayExperiment::colData(measurements) <- MultiAssayExperiment::colData(measurements)[colnames(MultiAssayExperiment::colData(measurements)) %in% sampleInfoTargets[, 2]]
  measurements <- wideFormat(measurements, colDataCols = seq_along(MultiAssayExperiment::colData(measurements)), check.names = FALSE, collapse = ':')
  measurements <- measurements[, -1, drop = FALSE] # Remove sample IDs.
  S4Vectors::mcols(measurements)[, "sourceName"] <- gsub("colDataCols", "clinical", S4Vectors::mcols(measurements)[, "sourceName"])
  colnames(S4Vectors::mcols(measurements))[1] <- "dataset"
  S4Vectors::mcols(measurements)[, "feature"] <- S4Vectors::mcols(measurements)[, "rowname"]
  missingIndices <- is.na(S4Vectors::mcols(measurements)[, "feature"])
  S4Vectors::mcols(measurements)[missingIndices, "feature"] <- colnames(measurements)[missingIndices]
  S4Vectors::mcols(measurements) <- S4Vectors::mcols(measurements)[, c("dataset", "feature")]

  plotFeatureClasses(measurements, classes, S4Vectors::mcols(measurements), groupBy, groupingName, showDatasetName = showDatasetName, ...)
})