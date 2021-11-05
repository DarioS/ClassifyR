setGeneric("selectionPlot", function(results, ...)
standardGeneric("selectionPlot"))

setMethod("selectionPlot", "list", 
          function(results,
                   comparison = "within", referenceLevel = NULL,
                   characteristicsList = list(x = "Classifier Name"), coloursList = list(), orderingList = list(), binsList = list(),
                   yMax = 100, fontSizes = c(24, 16, 12, 16), title = if(comparison[1] == "within") "Feature Selection Stability" else if(comparison == "size") "Feature Selection Size" else "Feature Selection Commonality",
                   yLabel = if(is.null(referenceLevel) && comparison != "size") "Common Features (%)" else if(comparison == "size") "Set Size" else paste("Common Features with", referenceLevel, "(%)"),
                   margin = grid::unit(c(1, 1, 1, 1), "lines"), rotate90 = FALSE, showLegend = TRUE, plot = TRUE, parallelParams = bpparam())
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")             
  if(!requireNamespace("scales", quietly = TRUE))
    stop("The package 'scales' could not be found. Please install it.")
  if(comparison == "within" && !is.null(referenceLevel))
    stop("'comparison' should not be \"within\" if 'referenceLevel' is not NULL.")              

  ggplot2::theme_set(ggplot2::theme_classic() + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA)))            
  
  allFeaturesList <- lapply(results, function(result)
  {
    .getFeaturesStrings(features(result))
  })

  if(!is.null(characteristicsList[['x']]))
    xValues <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == characteristicsList[['x']], "value"])
  if(!is.null(characteristicsList[["row"]]))
    rowValues <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == characteristicsList[["row"]], "value"])
  if(!is.null(characteristicsList[["column"]]))
    columnValues <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == characteristicsList[["column"]], "value"])
  if(comparison != "within")
    referenceVar <- sapply(results, function(result) result@characteristics[result@characteristics[, "characteristic"] == comparison, "value"])
  
  allCharacteristics <- c(unlist(characteristicsList), comparison)
  allCharacteristics <- setdiff(allCharacteristics, c("within", "size"))

  if(!is.null(referenceLevel) && !(referenceLevel %in% referenceVar))
    stop("Reference level is neither a level of the comparison factor nor is it NULL.")
  
  if(comparison == "within")
  {
    plotData <- do.call(rbind, bpmapply(function(result, featuresList)
    {
      percentOverlaps <- unlist(mapply(function(features, index)
      {
        otherFeatures <- featuresList[(index + 1):length(featuresList)]
        sapply(otherFeatures, function(other)
        {
          length(intersect(features, other)) / length(union(features, other)) * 100
        })
      }, featuresList[1:(length(featuresList) - 1)], 1:(length(featuresList) - 1), SIMPLIFY = FALSE))      

      characteristicsOrder <- match(allCharacteristics, result@characteristics[["characteristic"]])
      characteristicsList <- as.list(result@characteristics[["value"]])[characteristicsOrder]
      summaryTable <- data.frame(characteristicsList, overlap = percentOverlaps, check.names = FALSE)
      colnames(summaryTable)[1:length(allCharacteristics)] <- allCharacteristics
      summaryTable
    }, results, allFeaturesList, BPPARAM = parallelParams, SIMPLIFY = FALSE))
  } else if(comparison == "size")
  {
    plotData <- do.call(rbind, mapply(function(result, featuresList)
                {
                  setSizes <- sapply(featuresList, length)
                  characteristicsOrder <- match(allCharacteristics, result@characteristics[["characteristic"]])
                  characteristicsList <- as.list(result@characteristics[["value"]])[characteristicsOrder]
                  summaryTable <- data.frame(characteristicsList, size = setSizes, check.names = FALSE)
                  colnames(summaryTable)[1:length(allCharacteristics)] <- allCharacteristics
                  summaryTable
                }, results, allFeaturesList, SIMPLIFY = FALSE))

    plotData[, "size"] <- cut(plotData[, "size"], breaks = binsList[["setSizes"]])
    selectionSizes <- as.data.frame(table(plotData)) # Spaces will be destroyed.
    colnames(selectionSizes)[1:length(characteristicsList)] <- unlist(characteristicsList)
    groupingIDs <- apply(selectionSizes[, -match(c("size", "Freq"), colnames(selectionSizes)), drop = FALSE], 1, paste, collapse = '.')
    plotData <- do.call(rbind, by(selectionSizes, groupingIDs, function(dataSubset) {
                                  dataSubset[, "Freq"] <- dataSubset[, "Freq"] / sum(dataSubset[, "Freq"]) * 100
                                  dataSubset                                                   
                                  })
                        )
    
    plotData[, "Freq"] <- cut(plotData[, "Freq"], breaks = binsList[["frequencies"]])
    if(min(binsList[["frequencies"]]) == 0)
    {
      levelsRestore <- levels(plotData[, "Freq"])
      plotData[, "Freq"] <- as.character(plotData[, "Freq"])
      plotData[is.na(plotData[, "Freq"]), "Freq"] <- 0
      plotData[, "Freq"] <- factor(plotData[, "Freq"], levels = c(0, levelsRestore))
    }
  } else { # Commonality analysis.
    groupingFactor <- paste(if('x' %in% names(characteristicsList) && characteristicsList[['x']] != comparison) xValues,
                            if("row" %in% names(characteristicsList) && characteristicsList[["row"]] != comparison) rowValues,
                            if("column" %in% names(characteristicsList) && characteristicsList[["column"]] != comparison) columnValues, sep = " ")
    if(length(groupingFactor) == 0) groupingFactor <- rep("None", length(results))
    compareIndices <- split(1:length(results), groupingFactor)

    plotData <- do.call(rbind, bplapply(compareIndices, function(indicesSet)
    {
      if(is.null(referenceLevel))
      {
        indiciesCombinations <- lapply(1:length(indicesSet),
                                       function(index) c(indicesSet[index], indicesSet[-index]))
      } else { # Compare each factor level other than the reference level to the reference level.
        indiciesCombinations <- list(c(indicesSet[match(referenceLevel, referenceVar[indicesSet])],
                                       indicesSet[setdiff(1:length(indicesSet), match(referenceLevel, referenceVar[indicesSet]))]))
      }

      do.call(rbind, unname(lapply(indiciesCombinations, function(indiciesCombination)
      {
        aDataset <- results[[indiciesCombination[1]]]
        featuresList <- allFeaturesList[[indiciesCombination[1]]]
        otherResultIndices <- indiciesCombination[-1]
        otherResults <- results[otherResultIndices]
        
        overlapToOther <- lapply(otherResultIndices, function(otherIndex) # Other data sets to compare to.
        {
          unlist(lapply(featuresList, function(features) # List of features of a data set.
          {
            otherFeaturesList <- allFeaturesList[[otherIndex]]
            sapply(otherFeaturesList, function(otherFeatures) # List of features of another data set.
            {
              length(intersect(features, otherFeatures)) / length(union(features, otherFeatures)) * 100
            })
          }))
        })
        
        overlapToOther <- unlist(overlapToOther) # Convert all overlaps to a vector.
        if(is.null(referenceLevel))
        {
          characteristicsOrder <- match(allCharacteristics, aDataset@characteristics[["characteristic"]])
          characteristicsList <- as.list(aDataset@characteristics[["value"]])[characteristicsOrder]
          summaryTable <- data.frame(characteristicsList, overlap = overlapToOther)
          colnames(summaryTable)[1:length(allCharacteristics)] <- allCharacteristics
        } else { # Each other level has been compared to the reference level of the factor.
          otherSelections <- sapply(otherResults, length)
          
          summaryTable <- do.call(rbind, lapply(otherResults, function(otherResult)
          {
            selectTimes <- length(otherResults)
            do.call(cbind, lapply(allCharacteristics, function(characteristic)
              rep(otherResult@characteristics[otherResult@characteristics == characteristic, "value"], selectTimes)))
          }))
          colnames(summaryTable)[1:length(allCharacteristics)] <- allCharacteristics
        }      
        
        data.frame(summaryTable, overlap = overlapToOther, check.names = FALSE)
      })))
    }, BPPARAM = parallelParams))
  }
  rownames(plotData) <- NULL # Easier for viewing during maintenance.
  if(!"fillColours" %in% names(coloursList) && "fillColour" %in% names(characteristicsList))
  {
    if(characteristicsList[["fillColour"]] == "size")
    { # Automatically grey for zero.
      colours <- character()
      if(any(plotData[, "Freq"] == '0'))
        colours <- c("grey", scales::hue_pal()(length(unique(plotData[, "Freq"])) - 1))
      else
        colours <- scales::hue_pal()(length(unique(plotData[, "Freq"])))
      coloursList[["fillColours"]] <- colours
    } else {
      coloursList[["fillColours"]] <- scales::hue_pal()(length(unique(plotData[, characteristicsList[["fillColour"]]])))
    }
  }
  if(!"lineColours" %in% names(coloursList) && "lineColour" %in% names(characteristicsList))
    coloursList[["lineColours"]] <- scales::hue_pal(direction = -1)(length(unique(plotData[, characteristicsList[["lineColour"]]])))
  
  xLabel <- characteristicsList[['x']]
  xData <- plotData[, xLabel]
  if(length(orderingList) > 0) plotData <- .addUserLevels(plotData, orderingList)
  if(rotate90 == TRUE) plotData[, xLabel] <- factor(plotData[, xLabel], levels = rev(levels(plotData[, xLabel])))

  if(length(orderingList) > 0) plotData <- .addUserLevels(plotData, orderingList)
  characteristicsList <- lapply(characteristicsList, rlang::sym)
  
  if(comparison != "size")
  {
    legendPosition <- ifelse(showLegend == TRUE, "right", "none")
    selectionPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = !!characteristicsList[['x']], y = overlap, fill = !!characteristicsList[["fillColour"]], colour = !!characteristicsList[["lineColour"]])) +
                            ggplot2::coord_cartesian(ylim = c(0, yMax)) + ggplot2::xlab(xLabel) + ggplot2::ylab(yLabel) +
                            ggplot2::ggtitle(title) + ggplot2::theme(legend.position = legendPosition, axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5), plot.margin = margin)
    if(max(table(xData)) == 1) selectionPlot <- selectionPlot + ggplot2::geom_bar(stat = "identity") else selectionPlot <- selectionPlot + ggplot2::geom_violin()
  } else {
    selectionPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = !!characteristicsList[['x']], y = size)) +
                     ggplot2::geom_tile(ggplot2::aes(fill = Freq)) + ggplot2::ggtitle(title) + ggplot2::labs(x = xLabel, y = yLabel) + ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme(axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), plot.title = ggplot2::element_text(size = fontSizes[1], hjust = 0.5)) + ggplot2::guides(fill = ggplot2::guide_legend(title = "Frequency (%)"))
  }
  
  if(length(coloursList[["fillColours"]]) > 0)
    selectionPlot <- selectionPlot + ggplot2::scale_fill_manual(values = coloursList[["fillColours"]])
  if(length(coloursList[["lineColours"]]) > 0)
    selectionPlot <- selectionPlot + ggplot2::scale_colour_manual(values = coloursList[["lineColours"]])
  
  if(rotate90 == TRUE)
    selectionPlot <- selectionPlot + ggplot2::coord_flip(ylim = c(0, yMax))
  
  selectionPlot <- selectionPlot + ggplot2::facet_grid(ggplot2::vars(!!characteristicsList[["row"]]), ggplot2::vars(!!characteristicsList[["column"]])) + ggplot2::theme(strip.text = ggplot2::element_text(size = fontSizes[4]))

  if(plot == TRUE)
    print(selectionPlot)
  
  selectionPlot
})