setGeneric("rankingPlot", function(results, ...)
standardGeneric("rankingPlot"))

setMethod("rankingPlot", "list",
          function(results, topRanked = seq(10, 100, 10),
                   comparison = "within", referenceLevel = NULL,
                   characteristicsList = list(), orderingList = list(),
                   sizesList = list(lineWidth = 1, pointSize = 2, legendLinesPointsSize = 1,
                                    fonts = c(24, 16, 12, 12, 12, 16)),
                   lineColours = NULL, xLabelPositions = seq(10, 100, 10), yMax = 100,
                   title = if(comparison[1] == "within") "Feature Ranking Stability" else "Feature Ranking Commonality",
                   yLabel = if(is.null(referenceLevel)) "Average Common Features (%)" else paste("Average Common Features with", referenceLevel, "(%)"),
                   margin = grid::unit(c(1, 1, 1, 1), "lines"),
                   showLegend = TRUE, plot = TRUE, parallelParams = bpparam())
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")
  if(comparison == "within" && !is.null(referenceLevel))
    stop("'comparison' should not be \"within\" if 'referenceLevel' is not NULL.")
  nFeatures <- results[[1]]@selectResult@totalFeatures
  error <- character()
  if(max(topRanked) > nFeatures)
    error <- paste("'topRanked' is as high as", max(topRanked))
  if(max(xLabelPositions) > nFeatures)
    error <- paste(error, if(nchar(error) > 0) "and", "'xLabelPositions' is as high as", max(xLabelPositions))
  if(length(error) > 0)
  {
    error <- paste(error, "but there are only", nFeatures, "features in the data set.")
    stop(error)
  }
  ggplot2::theme_set(ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)))

  if(!is.null(characteristicsList[["lineColour"]]))
    lineColourValues <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == characteristicsList[["lineColour"]], "value"])
  if(!is.null(characteristicsList[["pointType"]]))
    pointTypeValues <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == characteristicsList[["pointType"]], "value"])
  if(!is.null(characteristicsList[["row"]]))
    rowValues <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == characteristicsList, "value"])
  if(!is.null(characteristicsList[["column"]]))
    columnValues <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == characteristicsList, "value"])
  if(comparison != "within")
    referenceVar <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == comparison, "value"])
  
  allCharacteristics <- c(unlist(characteristicsList), comparison)
  allCharacteristics <- setdiff(allCharacteristics, "within")
  
  if(!is.null(referenceLevel) && comparison != "within" && !(referenceLevel %in% referenceVar))
    stop("Reference level ", referenceLevel, " is neither a level of the comparison factor (",
          paste(referenceVar, collapse = ", "), ") nor is it NULL.")
  
  allRankedList <- lapply(results, function(result)
  {
    .getFeaturesStrings(result@selectResult@rankedFeatures)
  })
  
  if(comparison == "within")
  {
    plotData <- do.call(rbind, bpmapply(function(result, rankedList)
    {
      averageOverlap <- rowMeans(do.call(cbind, mapply(function(features, index)
      {
        otherFeatures <- rankedList[(index + 1):length(rankedList)]
        sapply(otherFeatures, function(other)
        {
          sapply(topRanked, function(top)
          {
            length(intersect(features[1:top], other[1:top])) / top * 100
          })
        })
      }, rankedList[1:(length(rankedList) - 1)], 1:(length(rankedList) - 1), SIMPLIFY = FALSE)))

      if(length(allCharacteristics) > 0)
      {
        characteristicsOrder <- match(allCharacteristics, result@characteristics[["characteristic"]])
        characteristicsList <- as.list(result@characteristics[["value"]])[characteristicsOrder]
        summaryTable <- data.frame(characteristicsList, top = topRanked, overlap = averageOverlap, check.names = FALSE)
        colnames(summaryTable)[1:length(allCharacteristics)] <- allCharacteristics
      } else {
        summaryTable <- data.frame(top = topRanked, overlap = averageOverlap)
      }
      summaryTable
    }, results, allRankedList, BPPARAM = parallelParams, SIMPLIFY = FALSE))
  } else { # Commonality analysis.
    groupingFactor <- paste(if("lineColour" %in% names(characteristicsList) && characteristicsList[["lineColour"]] != comparison) lineColourValues,
                            if("pointType" %in% names(characteristicsList) && characteristicsList[["pointType"]] != comparison) pointTypeValues,
                            if("row" %in% names(characteristicsList) && characteristicsList[["row"]] != comparison) rowValues,
                            if("column" %in% names(characteristicsList) && characteristicsList[["column"]] != comparison) columnValues, sep = " ")
    if(length(groupingFactor) == 0) groupingFactor <- rep("None", length(results))
    compareIndices <- split(1:length(results), groupingFactor)
    plotData <- do.call(rbind, unname(unlist(bplapply(compareIndices, function(indicesSet)
    {
      if(is.null(referenceLevel))
      {
        indiciesCombinations <- lapply(1:length(indicesSet),
                                       function(index) c(indicesSet[index], indicesSet[-index]))
      } else { # Compare each factor level other than the reference level to the reference level.
        indiciesCombinations <- list(c(indicesSet[match(referenceLevel, referenceVar[indicesSet])],
                                       indicesSet[setdiff(1:length(indicesSet), match(referenceLevel, referenceVar[indicesSet]))]))
      }

      unname(lapply(indiciesCombinations, function(indiciesCombination)
      {
        aDataset <- results[[indiciesCombination[1]]]
        rankedList <- allRankedList[[indiciesCombination[1]]]
        otherDatasetIndices <- indiciesCombination[-1]
        otherDatasets <- results[otherDatasetIndices]
        
        overlapToOther <- do.call(cbind, unname(lapply(otherDatasetIndices, function(otherIndex) # Other data sets to compare to.
        { 
          rowMeans(do.call(cbind, lapply(rankedList, function(rankings) # List of ranked features of a data set.
          {
            otherRankedList <- allRankedList[[otherIndex]]
            sapply(otherRankedList, function(otherRanked) # List of ranked features of another data set.
            {
              sapply(topRanked, function(top)
              {
                length(intersect(rankings[1:top], otherRanked[1:top])) / top * 100
              })          
            })
          })))
        })))

        if(is.null(referenceLevel))
        {
          overlapToOther <- rowMeans(overlapToOther)
          
          characteristicsOrder <- match(allCharacteristics, aDataset@characteristics[["characteristic"]])
          characteristicsList <- as.list(aDataset@characteristics[["value"]])[characteristicsOrder]
          summaryTable <- data.frame(characteristicsList, top = topRanked, overlap = overlapToOther, check.names = FALSE)
          colnames(summaryTable)[1:length(allCharacteristics)] <- allCharacteristics
          summaryTable
        } else { # Each other level has been compared to the reference level of the factor.
          overlapToOther <- as.numeric(overlapToOther) # Convert matrix of overlaps to a vector.
          topRankedAll <- rep(topRanked, length.out = length(overlapToOther))
          
          if(length(allCharacteristics) > 0)
            {
            summaryTable <- do.call(rbind, lapply(otherDatasets, function(otherDataset)
            {
              characteristicsOrder <- match(allCharacteristics, otherDataset@characteristics[["characteristic"]])
              characteristicsList <- as.list(otherDataset@characteristics[["value"]])[characteristicsOrder]
            }))
            colnames(summaryTable)[1:length(allCharacteristics)] <- allCharacteristics
            summaryTable <- summaryTable[rep(1:nrow(summaryTable), each = length(topRanked)), , drop = FALSE]
            data.frame(summaryTable, top = topRankedAll, overlap = overlapToOther,
                       check.names = FALSE)
          } else {
            data.frame(top = topRankedAll, overlap = overlapToOther)
          }
        }
      }))
    }, BPPARAM = parallelParams), recursive = FALSE)))
  }
  rownames(plotData) <- NULL # Easier for viewing during maintenance.
  if(is.null(lineColours) && "lineColour" %in% names(characteristicsList))
    lineColours <- scales::hue_pal()(length(unique(plotData[, characteristicsList[["lineColour"]]])))
  legendPosition <- ifelse(showLegend == TRUE, "right", "none")
  
  if(length(orderingList) > 0) plotData <- .addUserLevels(plotData, orderingList)
  if(length(characteristicsList) > 0) characteristicsList <- lapply(characteristicsList, rlang::sym)
  
  overlapPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = top, y = overlap, colour = !!characteristicsList[["lineColour"]], shape = !!characteristicsList[["pointType"]])) +
                          ggplot2::geom_line(size = sizesList[["lineWidth"]]) + ggplot2::geom_point(size = sizesList[["pointSize"]]) + ggplot2::scale_x_continuous(breaks = xLabelPositions, limits = range(xLabelPositions)) + ggplot2::coord_cartesian(ylim = c(0, yMax)) +
                          ggplot2::xlab("Top Features") + ggplot2::ylab(yLabel) + ggplot2::ggtitle(title) + ggplot2::scale_colour_manual(values = lineColours) +
                          ggplot2::theme(axis.title = ggplot2::element_text(size = sizesList[["fontSizes"]][2]), axis.text = ggplot2::element_text(colour = "black", size = sizesList[["fontSizes"]][3]), legend.position = legendPosition, legend.title = ggplot2::element_text(size = sizesList[["fontSizes"]][4]), legend.text = ggplot2::element_text(size = sizesList[["fontSizes"]][5]), plot.title = ggplot2::element_text(size = sizesList[["fontSizes"]][1], hjust = 0.5), plot.margin = margin) +
                          ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = sizesList[["legendLinesPointsSize"]])),
                                          shape = ggplot2::guide_legend(override.aes = list(size = sizesList[["legendLinesPointsSize"]])))
  
  overlapPlot <- overlapPlot + ggplot2::facet_grid(ggplot2::vars(!!characteristicsList[["row"]]), ggplot2::vars(!!characteristicsList[["column"]])) + ggplot2::theme(strip.text = ggplot2::element_text(size = sizesList[["fontSizes"]][6]))
  
  if(plot == TRUE)
    print(overlapPlot)
  
  overlapPlot
})
