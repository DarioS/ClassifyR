#' Plot Performance Measures for Various Classifications
#' 
#' Draws a graphical summary of a particular performance measure for a list of
#' classifications
#' 
#' If there are multiple values for a performance measure in a single result
#' object, it is plotted as a violin plot, unless \code{aggregate} is
#' \code{TRUE}, in which case the all predictions in a single result object are
#' considered simultaneously, so that only one performance number is
#' calculated, and a barchart is plotted.
#' 
#' @aliases performancePlot performancePlot,list-method
#' @param results A list of \code{\link{ClassifyResult}} objects.
#' @param aggregate A character vector of the levels of
#' \code{characteristicsList['x']} to aggregate to a single number by taking
#' the mean. This is particularly meaningful when the cross-validation is
#' leave-k-out, when k is small.
#' @param performanceName Default: "Balanced Accuracy". The name of the
#' performance measure to make comparisons of. This is one of the names printed
#' in the Performance Measures field when a \code{\link{ClassifyResult}} object is
#' printed, or if none are stored, the performance metric will be calculated.
#' @param characteristicsList A named list of characteristics. Each element's
#' name must be one of \code{"x"}, \code{"row"}, \code{"column"},
#' \code{fillColour}, or \code{fillLine}. The value of each element must be a
#' characteristic name, as stored in the \code{"characteristic"} column of the
#' results' characteristics table. Only \code{"x"} is mandatory.
#' @param coloursList A named list of plot aspects and colours for the aspects.
#' No elements are mandatory. If specified, each list element's name must be
#' either \code{"fillColours"} or \code{"lineColours"}. If a characteristic is
#' associated to fill or line by \code{characteristicsList} but this list is
#' empty, a palette of colours will be automaticaly chosen.
#' @param orderingList An optional named list. Any of the variables specified
#' to \code{characteristicsList} can be the name of an element of this list and
#' the value of the element is the order in which the factors should be
#' presented in, in case alphabetical sorting is undesirable.
#' @param yLimits The minimum and maximum value of the performance metric to
#' plot.
#' @param densityStyle Default: "box". Either \code{"violin"} for violin plot or
#' \code{"box"} for box plot.
#' @param fontSizes A vector of length 4. The first number is the size of the
#' title.  The second number is the size of the axes titles. The third number
#' is the size of the axes values. The fourth number is the font size of the
#' titles of grouped plots, if any are produced. In other words, when
#' \code{rowVariable} or \code{columnVariable} are not \code{NULL}.
#' @param title An overall title for the plot.
#' @param margin The margin to have around the plot.
#' @param rotate90 Logical. IF \code{TRUE}, the plot is horizontal.
#' @param showLegend If \code{TRUE}, a legend is plotted next to the plot. If
#' FALSE, it is hidden.
#' @param plot Logical. IF \code{TRUE}, a plot is produced on the current
#' graphics device.
#' @return An object of class \code{ggplot} and a plot on the current graphics
#' device, if \code{plot} is \code{TRUE}.
#' @author Dario Strbenac
#' @examples
#' 
#'   predicted <- data.frame(sample = sample(LETTERS[1:10], 80, replace = TRUE),
#'                           permutation = rep(1:2, each = 40),
#'                           class = factor(rep(c("Healthy", "Cancer"), 40)))
#'   actual <- factor(rep(c("Healthy", "Cancer"), each = 5))
#'   result1 <- ClassifyResult(DataFrame(characteristic = c("Data Set", "Selection Name", "Classifier Name",
#'                                                          "Cross-validation"),
#'                             value = c("Example", "t-test", "Differential Expression", "2 Permutations, 2 Folds")),
#'                             LETTERS[1:10], LETTERS[10:1], list(1:100, c(1:9, 11:101)),
#'                             list(c(1:3), c(2, 5, 6), 1:4, 5:8),
#'                             list(function(oracle){}), NULL, predicted, actual)
#'   result1 <- calcCVperformance(result1, "Macro F1")
#' 
#'   predicted <- data.frame(sample = sample(LETTERS[1:10], 80, replace = TRUE),
#'                           permutation = rep(1:2, each = 40),
#'                           class = factor(rep(c("Healthy", "Cancer"), 40)))
#'                                
#'   result2 <- ClassifyResult(DataFrame(characteristic = c("Data Set", "Selection Name", "Classifier Name",
#'                                                          "Cross-validation"),
#'                             value = c("Example", "Bartlett Test", "Differential Variability", "2 Permutations, 2 Folds")),
#'                             LETTERS[1:10], LETTERS[10:1], list(1:100, c(1:5, 11:105)),
#'                             list(c(1:3), c(4:6), c(1, 6, 7, 9), c(5:8)),
#'                             list(function(oracle){}), NULL, predicted, actual)
#'   result2 <- calcCVperformance(result2, "Macro F1")
#'   
#'   performancePlot(list(result1, result2), performanceName = "Macro F1",
#'                   title = "Comparison")
#' 
#' @importFrom rlang sym
#' @rdname performancePlot
#' @usage NULL
#' @export
setGeneric("performancePlot", function(results, ...) standardGeneric("performancePlot"))

#' @rdname performancePlot
#' @export
setMethod("performancePlot", "list", 
          function(results, performanceName = "Balanced Accuracy",
                   characteristicsList = list(x = "Classifier Name"), aggregate = character(), coloursList = list(), orderingList = list(),
                   densityStyle = c("box", "violin"), yLimits = NULL, fontSizes = c(24, 16, 12, 12), title = NULL,
                   margin = grid::unit(c(1, 1, 1, 1), "lines"), rotate90 = FALSE, showLegend = TRUE, plot = TRUE)
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")             
  if(!requireNamespace("scales", quietly = TRUE))
    stop("The package 'scales' could not be found. Please install it.")
  densityStyle <- match.arg(densityStyle)
  densityStyle <- ifelse(densityStyle == "box", ggplot2::geom_boxplot, ggplot2::geom_violin)
            
  ggplot2::theme_set(ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)))
  performanceNames <- unlist(lapply(results, function(result)
    if(!is.null(result@performance)) names(result@performance)))
  namesCounts <- table(performanceNames)
  commonNames <- names(namesCounts)[namesCounts == length(results)]
  if(!performanceName %in% commonNames)
  {
    warning(paste(performanceName, "not found in all elements of results. Calculating it now."))
    results <- lapply(results, function(result) calcCVperformance(result, performanceName))
  }
  
  if("dataset"%in%results[[1]]@characteristics$characteristic&!"dataset"%in%unlist(characteristicsList)){
    characteristicsList[["dataset"]] <- "dataset"
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
  
  plotData <- plotData[,!duplicated(colnames(plotData))]
  
  if("dataset"%in%colnames(plotData)&!"dataset"%in%unlist(characteristicsList[names(characteristicsList)!="dataset"])){
    if(length(unique(plotData$dataset))>1&!"fillColour"%in%names(characteristicsList)){
      characteristicsList[["fillColour"]] <- "dataset"
    }
    if(length(unique(plotData$dataset))>1&!"lineColour"%in%names(characteristicsList)&!"dataset"%in%unlist(characteristicsList[names(characteristicsList)!="dataset"])){
      characteristicsList[["lineColour"]] <- "dataset"
    }
  }  

  if("fillColour" %in% names(characteristicsList))
    if(!"fillColours" %in% names(coloursList)) coloursList[["fillColours"]] <- scales::hue_pal()(length(unique(plotData[, characteristicsList[["fillColour"]]])))
  if("lineColour" %in% names(characteristicsList))
    if(!"lineColours" %in% names(coloursList)) coloursList[["lineColours"]] <- scales::hue_pal(direction = -1)(length(unique(plotData[, characteristicsList[["lineColour"]]])))
  
  allCharacteristics <- unlist(characteristicsList)
  xLabel <- allCharacteristics['x']
  if(rotate90 == TRUE) plotData[, xLabel] <- factor(plotData[, xLabel], levels = rev(levels(plotData[, xLabel])))

  legendPosition <- ifelse(showLegend == TRUE, "right", "none")
  characteristicsList <- lapply(characteristicsList, rlang::sym)

  performancePlot <- ggplot2::ggplot() + 
                          ggplot2::ggtitle(title) + ggplot2::theme(legend.position = legendPosition, axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5), plot.margin = margin) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 2)

  if(!is.null(yLimits)) performancePlot <- performancePlot + ggplot2::coord_cartesian(ylim = yLimits)
  if("fillColour" %in% names(characteristicsList))
    performancePlot <- performancePlot + ggplot2::scale_fill_manual(values = coloursList[["fillColours"]])
  if("lineColour" %in% names(characteristicsList))
    performancePlot <- performancePlot + ggplot2::scale_colour_manual(values = coloursList[["lineColours"]])

  analysisGrouped <- split(plotData, plotData[, allCharacteristics])
  analysisGroupSizes <- sapply(analysisGrouped, nrow)
  if(any(analysisGroupSizes > 1))
  {
    multiPlotData <- do.call(rbind, analysisGrouped[analysisGroupSizes > 1])
    performancePlot <- performancePlot + densityStyle(data = multiPlotData, ggplot2::aes(x = !!characteristicsList[['x']], y = !!(rlang::sym(performanceName)), fill = !!characteristicsList[["fillColour"]], colour = !!characteristicsList[["lineColour"]]))
  }
  if(any(analysisGroupSizes == 1))
  {
    singlePlotData <- do.call(rbind, analysisGrouped[analysisGroupSizes == 1])
    performancePlot <- performancePlot + ggplot2::geom_bar(data = singlePlotData, stat = "identity", ggplot2::aes(x = !!characteristicsList[['x']], fill = !!characteristicsList[['fillColour']]),
                                                           colour = !!characteristicsList[['lineColour']])
  }
  
  if(!is.null(yLimits)) yLimits = c(0, 1)
  if(rotate90 == TRUE)
    performancePlot <- performancePlot + ggplot2::coord_flip(ylim = yLimits)
  
  performancePlot <- performancePlot + ggplot2::facet_grid(ggplot2::vars(!!characteristicsList[["row"]]), ggplot2::vars(!!characteristicsList[["column"]])) + ggplot2::theme(strip.text = ggplot2::element_text(size = fontSizes[4]))

  if(plot == TRUE)
    print(performancePlot)
  
  performancePlot
})
