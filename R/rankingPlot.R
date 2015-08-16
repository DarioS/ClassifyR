setGeneric("rankingPlot", function(results, ...)
{standardGeneric("rankingPlot")})

setMethod("rankingPlot", "list", 
          function(results, topRanked = seq(10, 100, 10),
                   comparison = c("within", "classificationName", "validation", "datasetName", "selectionName"),
                   referenceLevel = NULL,
                   lineColourVariable = c("validation", "datasetName", "classificationName", "selectionName", "None"),
                   lineColours = NULL, lineWidth = 1,
                   pointTypeVariable = c("datasetName", "classificationName", "validation", "selectionName", "None"),
                   pointSize = 2, legendLinesPointsSize = 1,
                   rowVariable = c("None", "datasetName", "classificationName", "validation", "selectionName"),
                   columnVariable = c("classificationName", "datasetName", "validation", "selectionName", "None"),
                   yMax = 100, fontSizes = c(24, 16, 12, 12, 12, 16),
                   title = if(comparison == "within") "Feature Ranking Stability" else "Feature Ranking Commonality",
                   xLabelPositions = seq(10, 100, 10),
                   yLabel = if(is.null(referenceLevel)) "Average Pairwise Common Features (%)" else paste("Average Common Features with", referenceLevel, "(%)"),
                   plot = TRUE, parallelParams = bpparam())
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")
  if(comparison == "within" && !is.null(referenceLevel))
    stop("'comparison' should not be \"within\" if 'referenceLevel' is not NULL.")            

  comparison <- match.arg(comparison)
  lineColourVariable <- match.arg(lineColourVariable)
  pointTypeVariable <- match.arg(pointTypeVariable)
  rowVariable <- match.arg(rowVariable)
  columnVariable <- match.arg(columnVariable)
  if(class(results[[1]]) == "ClassifyResult") resultsType <- "classification" else resultsType <- "selection"

  if(resultsType == "classification")
  {
    analyses <- sapply(results, function(result) result@classificationName)
    selections <- sapply(results, function(result) result@selectResult@selectionName)
    validations <- sapply(results, function(result) .validationText(result))
  } else { # Compare selections.
    analyses <- rep("Not classification", length(results))
    selections <- sapply(results, function(result) result@selectionName)
    validations <- rep("No cross-validation", length(results))
  }
  datasets <- sapply(results, function(result) result@datasetName)
  
  referenceVar <- switch(comparison, classificationName = analyses,
                         selectionName = selections,
                         datasetName = datasets,
                         validation = validations)
  if(!is.null(referenceLevel) && !(referenceLevel %in% referenceVar))
    stop("Reference level is neither a level of the comparison factor nor is it NULL.")
  
  if(resultsType == "classification")
  {
    results <- lapply(results, function(result)
    {
      if(result@validation[[1]] == "fold")  # Unlist the folds of Resample and Fold.
        result@selectResult@rankedFeatures <- unlist(result@selectResult@rankedFeatures, recursive = FALSE)
      result
    })
  }
  
  if(comparison == "within")
  {
    if(resultsType == "selection")
      stop("'comparison' should not be \"within\" for results that are not cross-validations.")
    plotData <- do.call(rbind, bplapply(results, function(result)
    {
      rankedFeatures <- result@selectResult@rankedFeatures
      averageOverlap <- rowMeans(do.call(cbind, mapply(function(features, index)
      {
        otherFeatures <- rankedFeatures[(index + 1):length(rankedFeatures)]
        sapply(otherFeatures, function(other)
        {
          sapply(topRanked, function(top)
          {
            length(intersect(features[1:top], other[1:top])) / top * 100
          })
        })
      }, rankedFeatures[1:(length(rankedFeatures) - 1)], 1:(length(rankedFeatures) - 1), SIMPLIFY = FALSE)))
      validationText <- .validationText(result)
      
      data.frame(dataset = rep(result@datasetName, length(topRanked)),
                 analysis = rep(result@classificationName, length(topRanked)),
                 selection = rep(result@selectResult@selectionName, length(topRanked)),
                 validation = rep(validationText, length(topRanked)),
                 top = topRanked,
                 overlap = averageOverlap)
    }, BPPARAM = parallelParams))
  } else { # Commonality analysis.
    compareIndices <- switch(comparison,
                             classificationName = split(1:length(results),
                                                        paste(validations, datasets, selections, sep = " ")),
                             validation = split(1:length(results), paste(analyses, datasets, selections, sep = " ")),
                             datasetName = split(1:length(results), paste(analyses, validations, selections, sep = " ")),
                             selectionName = split(1:length(results), paste(validations, datasets, sep = " "))
                       )
    
    plotData <- do.call(rbind, unname(unlist(bplapply(compareIndices, function(indicesSet)
    {
      if(is.null(referenceLevel))
      {
        indiciesCombinations <- lapply(1:length(indicesSet),
                                       function(index) c(indicesSet[index], indicesSet[-index]))
      } else { # Compare each factor level other than the reference level to the reference level.
        indiciesCombinations <- list(c(indicesSet[match(referenceLevel, referenceVar)],
                                       indicesSet[setdiff(length(indicesSet), match(referenceLevel, referenceVar))]))
      }

      unname(lapply(indiciesCombinations, function(indiciesCombination)
      {
        aDataset <- results[[indiciesCombination[1]]]
        otherDatasets <- results[indiciesCombination[-1]]
        if(resultsType == "classification")
          featuresList <- aDataset@selectResult@rankedFeatures
        else
          featuresList <- aDataset@rankedFeatures
        
        overlapToOther <- do.call(cbind, unname(lapply(otherDatasets, function(anotherDataset) # Other datasets to compare to.
        { 
          rowMeans(do.call(cbind, lapply(featuresList, function(features) # List of features of a dataset.
          {
            if(resultsType == "classification")
              otherFeaturesList <- anotherDataset@selectResult@rankedFeatures
            else
              otherFeaturesList <- anotherDataset@rankedFeatures

            sapply(otherFeaturesList, function(otherFeatures) # List of features of another dataset.
            {
              sapply(topRanked, function(top)
              {
                length(intersect(features[1:top], otherFeatures[1:top])) / top * 100
              })          
            })
          })))
        })))
        
        if(is.null(referenceLevel))
        {
          overlapToOther <- rowMeans(overlapToOther)
          datasetText <- rep(aDataset@datasetName, length(topRanked))
          
          if(resultsType == "classification")
          {
            selectionText <- rep(aDataset@selectResult@selectionName, length(topRanked))
            analysisText <- rep(aDataset@classificationName, length(topRanked))
            validationText <- .validationText(aDataset)
          } else { # For standalone feature selection, there is no classification.
            selectionText <- rep(aDataset@selectionName, length(topRanked))
            analysisText <- "No classification"
            validationText <- "No cross-validation"
          }
        } else { # Each other level has been compared to the reference level of the factor.
          overlapToOther <- as.numeric(overlapToOther) # Convert matrix of overlaps to a vector.
          datasetText <- rep(sapply(otherDatasets, function(dataset) dataset@datasetName), each = length(topRanked))
          selectionText <- rep(sapply(otherDatasets, function(dataset) dataset@selectionName),
                                                                       each = length(topRanked))
          if(resultsType == "classification")
          {
            analysisText <- rep(sapply(otherDatasets, function(dataset) dataset@classificationName),
                                each = length(topRanked))
            validationText <- rep(sapply(otherDatasets, function(dataset) .validationText(dataset)),
                                  each = length(topRanked))
          } else { # For standalone feature selection, there is no classification.
            analysisText <- "No classification"
            validationText <- "No cross-validation"
          }
        }
        topRankedAll <- rep(topRanked, length.out = length(overlapToOther))
        
        data.frame(dataset = datasetText,
                   analysis = analysisText,
                   selection = selectionText,
                   validation = validationText,
                   top = topRankedAll,
                   overlap = overlapToOther)
      }))
    }, BPPARAM = parallelParams), recursive = FALSE)))
  }

  # Order factors in which they appeared in the user's list.
  plotData[, "dataset"] <- factor(plotData[, "dataset"], levels = unique(datasets))
  plotData[, "analysis"] <- factor(plotData[, "analysis"], levels = unique(analyses))
  plotData[, "selection"] <- factor(plotData[, "selection"], levels = unique(selections))
  plotData[, "validation"] <- factor(plotData[, "validation"], levels = unique(validations))
  
  if(is.null(lineColours) && lineColourVariable != "None")
    lineColours <- scales::hue_pal()(switch(lineColourVariable, validation = length(unique(plotData[, "validation"])), datasetName = length(unique(plotData[, "dataset"])), classificationName = length(unique(plotData[, "analysis"])), selectionName = length(unique(plotData[, "selection"]))))

  overlapPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = top, y = overlap,
                          colour = switch(lineColourVariable, validation = validation, datasetName = dataset, classificationName = analysis, selectionName = selection, None = NULL),
                          shape = switch(pointTypeVariable, validation = validation, datasetName = dataset, classificationName = analysis, selectionName = selection, None = NULL)), environment = environment()) +
                          ggplot2::geom_line(size = lineWidth) + ggplot2::geom_point(size = pointSize) + ggplot2::scale_x_continuous(breaks = xLabelPositions, limits = range(xLabelPositions)) + ggplot2::scale_y_continuous(limits = c(0, yMax)) +
                          ggplot2::xlab("Top Features") + ggplot2::ylab(yLabel) +
                          ggplot2::ggtitle(title) + ggplot2::labs(colour = switch(lineColourVariable, validation = "Validation", datasetName = "Dataset", classificationName = "Analysis", classificationName = "Analysis", selectionName = "Feature\nSelection"), shape = switch(pointTypeVariable, validation = "Validation", datasetName = "Dataset", classificationName = "Analysis", selectionName = "Feature\nSelection")) +
                          ggplot2::theme(axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), legend.title = ggplot2::element_text(size = fontSizes[4]), legend.text = ggplot2::element_text(size = fontSizes[5]), plot.title = ggplot2::element_text(size = fontSizes[1]), plot.margin = grid::unit(c(0, 0, 1, 0), "lines"), legend.margin = grid::unit(-1, "lines")) +
                          ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = legendLinesPointsSize)),
                                          shape = ggplot2::guide_legend(override.aes = list(size = legendLinesPointsSize)))
  
  if(!is.null(lineColours))
    overlapPlot <- overlapPlot + ggplot2::scale_colour_manual(values = lineColours)
  
  if(rowVariable != "None" || columnVariable != "None")
    overlapPlot <- overlapPlot + ggplot2::facet_grid(paste(if(rowVariable != "None") switch(rowVariable, validation = "validation", datasetName = "dataset", classificationName = "analysis"), "~", if(columnVariable != "None") switch(columnVariable, validation = "validation", datasetName = "dataset", classificationName = "analysis"))) + ggplot2::theme(strip.text = ggplot2::element_text(size = fontSizes[6]))
  
  if(plot == TRUE)
    print(overlapPlot)
  
  overlapPlot
})
