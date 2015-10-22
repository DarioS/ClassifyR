setGeneric("selectionPlot", function(results, ...)
{standardGeneric("selectionPlot")})

setMethod("selectionPlot", "list", 
          function(results,
                   comparison = c("within", "size", "classificationName", "validation", "datasetName", "selectionName"),
                   referenceLevel = NULL,
                   xVariable = c("classificationName", "datasetName", "validation", "selectionName"),
                   boxFillColouring = c("classificationName", "size", "datasetName", "validation", "selectionName", "None"),
                   boxFillColours = NULL,
                   boxFillBinBoundaries = NULL,
                   setSizeBinBoundaries = NULL,
                   boxLineColouring = c("validation", "classificationName", "datasetName", "selectionName", "None"),
                   boxLineColours = NULL,
                   rowVariable = c("None", "validation", "datasetName", "classificationName", "selectionName"),
                   columnVariable = c("datasetName", "classificationName", "validation", "selectionName", "None"),
                   yMax = 100, fontSizes = c(24, 16, 12, 16), title = if(comparison[1] == "within") "Feature Selection Stability" else if(comparison == "size") "Feature Selection Size" else "Feature Selection Commonality",
                   xLabel = "Analysis", yLabel = if(is.null(referenceLevel) && comparison != "size") "Common Features (%)" else if(comparison == "size") "Set Size" else paste("Common Features with", referenceLevel, "(%)"),
                   margin = grid::unit(c(0, 1, 1, 0), "lines"), rotate90 = FALSE, plot = TRUE, parallelParams = bpparam())
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")             
  if(!requireNamespace("scales", quietly = TRUE))
    stop("The package 'scales' could not be found. Please install it.")
  if(comparison == "within" && !is.null(referenceLevel))
    stop("'comparison' should not be \"within\" if 'referenceLevel' is not NULL.")                

  comparison <- match.arg(comparison)
  xVariable <- match.arg(xVariable)
  boxFillColouring <- match.arg(boxFillColouring)
  boxLineColouring <- match.arg(boxLineColouring)
  rowVariable <- match.arg(rowVariable)
  columnVariable <- match.arg(columnVariable)
  if(class(results[[1]]) == "ClassifyResult") resultsType <- "classification" else resultsType <- "selection"
  
  if(resultsType == "classification")
  {
    results <- lapply(results, function(result)
    {
      if(result@validation[[1]] == "fold")  # Unlist the folds of Resample and Fold.
        result@selectResult@chosenFeatures <- unlist(features(result), recursive = FALSE)
      result
    })
  }

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
  
  if(comparison == "within")
  {
    if(resultsType == "selection")
      stop("'comparison' should not be \"within\" for results that are not cross-validations.")
    plotData <- do.call(rbind, bplapply(results, function(result)
    {
      chosenFeatures <- features(result)
      percentOverlaps <- unlist(mapply(function(features, index)
      {
        otherFeatures <- chosenFeatures[(index + 1):length(chosenFeatures)]
        sapply(otherFeatures, function(other)
        {
          length(intersect(features, other)) / length(union(features, other)) * 100
        })
      }, chosenFeatures[1:(length(chosenFeatures) - 1)], 1:(length(chosenFeatures) - 1), SIMPLIFY = FALSE))      
      validationText <- .validationText(result)
      
      data.frame(dataset = rep(result@datasetName, length(percentOverlaps)),
                 analysis = rep(result@classificationName, length(percentOverlaps)),
                 selection = rep(result@selectResult@selectionName, length(percentOverlaps)),
                 validation = rep(validationText, length(percentOverlaps)),
                 overlap = percentOverlaps)
    }, BPPARAM = parallelParams))
  } else if(comparison == "size")
  {
    plotData <- do.call(rbind, lapply(results, function(result)
                {
                  selectedFeatures <- features(result)
                  setSizes <- sapply(selectedFeatures, length)
                  validation <- .validationText(result)
                  data.frame(analysis = rep(result@classificationName, length(setSizes)),
                             dataset = rep(result@datasetName, length(setSizes)),
                             selection = rep(result@selectResult@selectionName, length(setSizes)),
                             validation = rep(validation, length(setSizes)),
                             size = setSizes)
                }))

    plotData[, "size"] <- cut(plotData[, "size"], breaks = setSizeBinBoundaries, include.lowest = TRUE)
    selectionSizes <- as.data.frame(table(plotData[, 1:5]))

    plotData <- do.call(rbind, by(selectionSizes, apply(selectionSizes[, c("analysis", "dataset", "selection", "validation")], 1, paste, collapse = '.'),
                                              function(dataSubset) {
                                                dataSubset[, "Freq"] <- dataSubset[, "Freq"] / sum(dataSubset[, "Freq"]) * 100
                                                dataSubset                                                   
                                              }))
    plotData[, "Freq"] <- cut(plotData[, "Freq"], breaks = boxFillBinBoundaries, include.lowest = TRUE)
  } else { # Commonality analysis.
    groupingVariablesValues <- setdiff(c(xVariable, boxFillColouring, boxLineColouring, rowVariable, columnVariable), comparison)
    groupingFactor <- paste(if("datasetName" %in% groupingVariablesValues) datasets,
                            if("classificationName" %in% groupingVariablesValues) analyses,
                            if("selectionName" %in% groupingVariablesValues) selections,
                            if("validation" %in% groupingVariablesValues) validations,
                            sep = " ")
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
        otherDatasets <- results[indiciesCombination[-1]]
        if(resultsType == "classification")
          featuresList <- features(aDataset)
        else
          featuresList <- aDataset@chosenFeatures        
        
        overlapToOther <- lapply(otherDatasets, function(anotherDataset) # Other datasets to compare to.
        {
          unlist(lapply(featuresList, function(features) # List of features of a dataset.
          {
            if(resultsType == "classification")
              otherFeaturesList <- features(anotherDataset)
            else
              otherFeaturesList <- anotherDataset@chosenFeatures
            
            sapply(otherFeaturesList, function(otherFeatures) # List of features of another dataset.
            {
              length(intersect(features, otherFeatures)) / length(union(features, otherFeatures)) * 100
            })
          }))
        })
        
        if(is.null(referenceLevel))
        {
          overlapToOther <- unlist(overlapToOther)
          datasetText <- rep(aDataset@datasetName, length(overlapToOther))
          
          if(resultsType == "classification")
          {
            selectionText <- rep(aDataset@selectResult@selectionName, length(overlapToOther))
            analysisText <- rep(aDataset@classificationName, length(overlapToOther))
            validationText <- .validationText(aDataset)
          } else { # For standalone feature selection, there is no classification.
            selectionText <- rep(aDataset@selectionName, length(overlapToOther))
            analysisText <- "No classification"
            validationText <- "No cross-validation"
          }
        } else { # Each other level has been compared to the reference level of the factor.
          otherSelections <- sapply(otherDatasets, length)
          datasetText <- rep(sapply(otherDatasets, function(dataset) dataset@datasetName), otherSelections)
          selectionText <- rep(sapply(otherDatasets, function(dataset) if(resultsType == "classification") dataset@selectResult@selectionName else dataset@selectionName), otherSelections)
          
          if(resultsType == "classification")
          {
            analysisText <- rep(sapply(otherDatasets, function(dataset) dataset@classificationName),
                                otherSelections)
            validationText <- rep(sapply(otherDatasets, function(dataset) .validationText(dataset)),
                                  otherSelections)
          } else { # For standalone feature selection, there is no classification.
            analysisText <- "No classification"
            validationText <- "No cross-validation"
          }
          overlapToOther <- as.numeric(overlapToOther) # Convert matrix of overlaps to a vector.
        }      
        
        data.frame(dataset = datasetText,
                   analysis = analysisText,
                   selection = selectionText,
                   validation = validationText,
                   overlap = overlapToOther)
      })))
    }, BPPARAM = parallelParams))
  }
  rownames(plotData) <- NULL # Easier for viewing during maintenance.
  
  if(is.null(boxFillColours) && boxFillColouring != "None")
    boxFillColours <- scales::hue_pal()(switch(boxFillColouring, validation = length(unique(plotData[, "validation"])), datasetName = length(unique(plotData[, "dataset"])), classificationName = length(unique(plotData[, "analysis"])), selectionName = length(unique(plotData[, "selection"])), size = length(unique(plotData[, "Freq"]))))
  if(is.null(boxLineColours) && boxLineColouring != "None")
    boxLineColours <- scales::hue_pal(direction = -1)(switch(boxLineColouring, validation = length(unique(plotData[, "validation"])), datasetName = length(unique(plotData[, "dataset"])), classificationName = length(unique(plotData[, "analysis"])), selectionName = length(unique(plotData[, "selection"]))))
  
  # Order factors in which they appeared in the user's list.
  plotData[, "dataset"] <- factor(plotData[, "dataset"], levels = unique(datasets))
  plotData[, "analysis"] <- factor(plotData[, "analysis"], levels = unique(analyses))
  plotData[, "selection"] <- factor(plotData[, "selection"], levels = unique(selections))
  plotData[, "validation"] <- factor(plotData[, "validation"], levels = unique(validations))

  if(rotate90 == TRUE)
  {
    switch(xVariable, 
           validation = plotData[, "validation"] <- factor(plotData[, "validation"], levels = rev(levels(plotData[, "validation"]))),
           datasetName = plotData[, "dataset"] <- factor(plotData[, "dataset"], levels = rev(levels(plotData[, "dataset"]))),
           classificationName = plotData[, "analysis"] <- factor(plotData[, "analysis"], levels = rev(levels(plotData[, "analysis"]))))
           selection = plotData[, "selection"] <- factor(plotData[, "selection"], levels = rev(levels(plotData[, "selection"])))
  }
  
  if(comparison != "size")
  {
    selectionPlot <- ggplot2::ggplot(plotData, ggplot2::aes_string(x = switch(xVariable, validation = "validation", datasetName = "dataset", classificationName = "analysis", selectionName = "selection"), y = "overlap",
                            fill = switch(boxFillColouring, validation = "validation", datasetName = "dataset", classificationName = "analysis", selectionName = "selection", None = NULL), colour = switch(boxLineColouring, validation = "validation", datasetName = "dataset", classificationName = "analysis", selectionName = "selection", None = NULL)), environment = environment()) +
                            ggplot2::scale_y_continuous(limits = c(0, yMax)) + ggplot2::xlab(xLabel) + ggplot2::ylab(yLabel) +
                            ggplot2::ggtitle(title) + ggplot2::theme(legend.position = "none", axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), plot.title = ggplot2::element_text(size = fontSizes[1]), plot.margin = margin)
    
    xData <- switch(xVariable, validation = plotData[, "validation"], datasetName = plotData[, "dataset"], classificationName = plotData[, "analysis"], selectionName = plotData[, "selection"])
    if(max(table(xData)) == 1) selectionPlot <- selectionPlot + ggplot2::geom_bar(stat = "identity") else selectionPlot <- selectionPlot + ggplot2::geom_boxplot()
  } else {
    selectionPlot <- ggplot2::ggplot(plotData, ggplot2::aes_string(x = switch(xVariable, validation = "validation", datasetName = "dataset", classificationName = "analysis", selectionName = "selection"), y = "size"), environment = environment()) +
                     ggplot2::geom_tile(ggplot2::aes(fill = Freq)) + ggplot2::ggtitle(title) + ggplot2::labs(x = xLabel, y = yLabel) + ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme(axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), plot.title = ggplot2::element_text(size = fontSizes[1])) + ggplot2::guides(fill = ggplot2::guide_legend(title = "Frequency (%)"))
  }

  if(!is.null(boxFillColours))
    selectionPlot <- selectionPlot + ggplot2::scale_fill_manual(values = boxFillColours)
  if(!is.null(boxLineColours))
    selectionPlot <- selectionPlot + ggplot2::scale_colour_manual(values = boxLineColours)
  
  if(rotate90 == TRUE)
    selectionPlot <- selectionPlot + ggplot2::coord_flip()
  
  if(rowVariable != "None" || columnVariable != "None")
    selectionPlot <- selectionPlot + ggplot2::facet_grid(paste(if(rowVariable != "None") switch(rowVariable, validation = "validation", datasetName = "dataset", classificationName = "analysis", selectionName = "selection"), "~", if(columnVariable != "None") switch(columnVariable, validation = "validation", datasetName = "dataset", classificationName = "analysis", selectionName = "selection"))) + ggplot2::theme(strip.text = ggplot2::element_text(size = fontSizes[4]))

  if(plot == TRUE)
    print(selectionPlot)
  
  selectionPlot
})