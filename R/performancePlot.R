setGeneric("performancePlot", function(results, ...)
standardGeneric("performancePlot"))

setMethod("performancePlot", "list", 
          function(results, performanceName = NULL,
                   characteristicsList = list(x = "Classifier Name"), aggregate = character(), coloursList = list(), orderingList = list(),
                   yLimits = c(0, 1), fontSizes = c(24, 16, 12, 12), title = NULL,
                   margin = grid::unit(c(1, 1, 1, 1), "lines"), rotate90 = FALSE, showLegend = TRUE, plot = TRUE)
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")             
  if(!requireNamespace("scales", quietly = TRUE))
    stop("The package 'scales' could not be found. Please install it.")
  
  ggplot2::theme_set(ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)))
  
  performanceNames <- lapply(results, function(result)
    if(!is.null(result@performance)) names(result@performance))
  namesCounts <- table(performanceNames)
  commonNames <- names(namesCounts)[namesCounts == length(results)]
  if(is.null(performanceName))
  {
    stop("Please specify a performance measure to plot.",
         "Calculated ones for all results are: ", paste(commonNames, collapse = ", "))
  }
  
  plotData <- do.call(rbind, mapply(function(result, index)
                    {
                      if(!performanceName %in% names(result@performance))
                        stop(performanceName, " not calculated for element ", index, " of results list.")
                      row <- result@characteristics[, "characteristic"] == characteristicsList[["x"]] 
                      if(any(row) && result@characteristics[row, "value"] %in% aggregate)
                        performance <- mean(result@performance[[performanceName]])
                      else
                        performance <- result@performance[[performanceName]]
                      rows <- match(unlist(characteristicsList), result@characteristics[, "characteristic"])
                      summaryTable <- data.frame(as.list(result@characteristics[rows, "value"]), performance)
                      colnames(summaryTable) <- c(characteristicsList, performanceName)
                      summaryTable
                    }, results, 1:length(results), SIMPLIFY = FALSE))

  if("fillColour" %in% names(characteristicsList))
    if(!"fillColours" %in% names(coloursList)) coloursList[["fillColours"]] <- scales::hue_pal()(length(unique(plotData[, xLabel])))
  if("lineColour" %in% names(characteristicsList))
    if(!"lineColours" %in% names(coloursList)) coloursList[["lineColours"]] <- scales::hue_pal(direction = -1)(length(unique(plotData[, xLabel])))
  
  allCharacteristics <- unlist(characteristicsList)
  xLabel <- allCharacteristics['x']
  if(rotate90 == TRUE) plotData[, xLabel] <- factor(plotData[, xLabel], levels = rev(levels(plotData[, xLabel])))

  legendPosition <- ifelse(showLegend == TRUE, "right", "none")
  characteristicsList <- lapply(characteristicsList, rlang::sym)
  performancePlot <- ggplot2::ggplot() + ggplot2::coord_cartesian(ylim = yLimits) +
                          ggplot2::ggtitle(title) + ggplot2::theme(legend.position = legendPosition, axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5), plot.margin = margin)

  if("fillColour" %in% names(characteristicsList))
    performancePlot <- performancePlot + ggplot2::scale_fill_manual(values = coloursList[["fillColours"]])
  if("lineColour" %in% names(characteristicsList))
    performancePlot <- performancePlot + ggplot2::scale_fill_manual(values = coloursList[["lineColours"]])

  analysisGrouped <- split(plotData, plotData[, allCharacteristics])
  analysisGroupSizes <- sapply(analysisGrouped, nrow)
  if(any(analysisGroupSizes > 1))
  {
    multiPlotData <- do.call(rbind, analysisGrouped[analysisGroupSizes > 1])
    performancePlot <- performancePlot + ggplot2::geom_violin(data = multiPlotData, ggplot2::aes(x = !!characteristicsList[['x']], y = !!(rlang::sym(performanceName)), fill = !!characteristicsList[["fillColour"]], colour = !!characteristicsList[["lineColour"]]))
  }
  if(any(analysisGroupSizes == 1))
  {
    singlePlotData <- do.call(rbind, analysisGrouped[analysisGroupSizes == 1])
    performancePlot <- performancePlot + ggplot2::geom_bar(data = singlePlotData, stat = "identity", ggplot2::aes(x = !!characteristicsList[['x']], fill = !!characteristicsList[['fillColour']]),
                                                           colour = !!characteristicsList[['lineColour']])
  }
  
  if(rotate90 == TRUE)
    performancePlot <- performancePlot + ggplot2::coord_flip(ylim = yLimits)
  
  performancePlot <- performancePlot + ggplot2::facet_grid(ggplot2::vars(!!characteristicsList[["row"]]), ggplot2::vars(!!characteristicsList[["column"]])) + ggplot2::theme(strip.text = ggplot2::element_text(size = fontSizes[4]))

  if(plot == TRUE)
    print(performancePlot)
  
  performancePlot
})