setGeneric("selectionPlot", function(results, ...)
{standardGeneric("selectionPlot")})

setMethod("selectionPlot", "list", 
          function(results,
                   comparison = c("within", "classificationName", "validation", "datasetName"),
                   xVariable = c("classificationName", "datasetName", "validation"),
                   boxFillColouring = c("classificationName", "datasetName", "validation", "None"),
                   boxFillColours = NULL,
                   boxLineColouring = c("validation", "classificationName", "validation", "None"),
                   boxLineColours = NULL,
                   rowVariable = c("None", "validation", "datasetName", "classificationName"),
                   columnVariable = c("datasetName", "classificationName", "validation", "None"),
                   yMax = 100, fontSizes = c(24, 16, 12), title = "Feature Selection Stability",
                   xLabel = "Analysis", yLabel = "Average Pairwise Common Features (%)",
                   margin = grid::unit(c(0, 1, 1, 0), "lines"), rotate90 = FALSE, plot = TRUE, parallelParams = bpparam())
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")             
  if(!requireNamespace("scales", quietly = TRUE))
    stop("The package 'scales' could not be found. Please install it.")             

  comparison <- match.arg(comparison)
  xVariable <- match.arg(xVariable)
  boxFillColouring <- match.arg(boxFillColouring)
  boxLineColouring <- match.arg(boxLineColouring)
  rowVariable <- match.arg(rowVariable)
  columnVariable <- match.arg(columnVariable)

  results <- lapply(results, function(result) # Unlist the folds of Resample and Fold.
  {
    if(is.list(result@chosenFeatures[[1]]))
      result@chosenFeatures <- unlist(result@chosenFeatures, recursive = FALSE)
    result
  })
  
  if(comparison == "within")
  {  
    plotData <- do.call(rbind, bplapply(results, function(result)
    {
  
      percentOverlaps <- unlist(mapply(function(features, index)
      {
        otherFeatures <- chosenFeatures[(index + 1):length(chosenFeatures)]
        sapply(otherFeatures, function(other)
        {
          length(intersect(features, other)) / max(length(features), length(other)) * 100
        })
      }, chosenFeatures[1:(length(chosenFeatures) - 1)], 1:(length(chosenFeatures) - 1), SIMPLIFY = FALSE))
      
      data.frame(dataset = rep(result@datasetName, length(percentOverlaps)),
                 analysis = rep(result@classificationName, length(percentOverlaps)),
                 validation = rep(result@validation[[1]], length(percentOverlaps)),
                 overlap = percentOverlaps)
    }, BPPARAM = parallelParams))
  } else { # Commonality analysis.
    analyses <- sapply(results, function(result) result@classificationName) 
    validations <- sapply(results, function(result) paste(unlist(result@validation), collapse = " "))
    datasets <- sapply(results, function(result) result@datasetName)
    
    if(comparison == "classificationName")
      compareIndices <- split(1:length(results), paste(validations, datasets, sep = " "))
    else if(comparison == "validation")
      compareIndices <- split(1:length(results), paste(analyses, datasets, sep = " "))
    else
      compareIndices <- split(1:length(results), paste(analyses, validations, sep = " "))
    
    plotData <- do.call(rbind, unname(unlist(lapply(compareIndices, function(indicesSet)
    {
      if(length(indicesSet) == 2)
        indiciesCombinations <- list(indicesSet)
      else
        indiciesCombinations <- lapply(1:length(indicesSet), function(index) c(indicesSet[index], indicesSet[-index]))      
      unname(lapply(indiciesCombinations, function(indiciesCombination)
      {
        aDataset <- results[[indiciesCombination[1]]]
        otherDatasets <- results[indiciesCombination[-1]]   
        allOverlaps <- unlist(lapply(aDataset@chosenFeatures, function(features)
        {
          sapply(otherDatasets, function(anotherDataset)
          {
            sapply(anotherDataset@chosenFeatures, function(otherFeatures)
            {
              length(intersect(features, otherFeatures)) / max(length(features), length(otherFeatures)) * 100
            })
          })
        }))

        data.frame(dataset = rep(aDataset@datasetName, length(allOverlaps)),
                   analysis = rep(aDataset@classificationName, length(allOverlaps)),
                   validation = rep(aDataset@validation[[1]], length(allOverlaps)),
                   overlap = allOverlaps)
      }))
    }), recursive = FALSE)))
  }

  if(boxFillColouring != "None")
    if(!is.null(boxFillColours)) boxFillColours else boxFillColours <- scales::hue_pal()(switch(boxFillColouring, validation = length(unique(plotData[, "validation"])), datasetName = length(unique(plotData[, "dataset"])), classificationName = length(unique(plotData[, "analysis"]))))
  if(boxLineColouring != "None")
    if(!is.null(boxLineColours)) boxLineColours else boxLineColours <- scales::hue_pal(direction = -1)(switch(boxLineColouring, validation = length(unique(plotData[, "validation"])), datasetName = length(unique(plotData[, "dataset"])), classificationName = length(unique(plotData[, "analysis"]))))
  
  validationOrder <- unique(sapply(results, function(result) result@validation[[1]]))
  validationLabels <- c("Resample and Fold", "Resample and Split", "Leave Out")
  newValidLevels <- sapply(validationOrder, function(validation) grep(validation, validationLabels, ignore.case = TRUE))
  
  plotData[, "dataset"] <- factor(plotData[, "dataset"], levels = unique(sapply(results, function(result) result@datasetName)))
  plotData[, "analysis"] <- factor(plotData[, "analysis"], levels = unique(sapply(results, function(result) result@classificationName)))
  plotData[, "validation"] <- gsub("fold", "Resample and Fold", plotData[, "validation"])
  plotData[, "validation"] <- gsub("split", "Resample and Split", plotData[, "validation"])
  plotData[, "validation"] <- gsub("leave", "Leave Out", plotData[, "validation"])
  plotData[, "validation"] <- factor(plotData[, "validation"], levels = validationLabels[newValidLevels])
  if(rotate90 == TRUE)
  {
    switch(xVariable, validation = plotData[, "validation"] <- factor(plotData[, "validation"], levels = rev(levels(plotData[, "validation"]))),
           datasetName = plotData[, "dataset"] <- factor(plotData[, "dataset"], levels = rev(levels(plotData[, "dataset"]))),
           classificationName = plotData[, "analysis"] <- factor(plotData[, "analysis"], levels = rev(levels(plotData[, "analysis"]))))
  }
  
  selectionPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = switch(xVariable, validation = validation, datasetName = dataset, classificationName = analysis), y = overlap,
                          fill = switch(boxFillColouring, validation = validation, datasetName = dataset, classificationName = analysis, None = NULL), colour = switch(boxLineColouring, validation = validation, datasetName = dataset, classificationName = analysis, None = NULL)), environment = environment()) +
                          ggplot2::geom_boxplot() + ggplot2::scale_y_continuous(limits = c(0, yMax)) + ggplot2::scale_fill_manual(values = boxFillColours) + ggplot2::scale_colour_manual(values = boxLineColours) + ggplot2::xlab(xLabel) + ggplot2::ylab(yLabel) +
                          ggplot2::ggtitle(title) + ggplot2::theme(legend.position = "none", axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), plot.title = ggplot2::element_text(size = fontSizes[1]), plot.margin = margin)
  if(rotate90 == TRUE)
    selectionPlot <- selectionPlot + ggplot2::coord_flip()
  
  if(rowVariable != "None" || columnVariable != "None")
    selectionPlot <- selectionPlot + ggplot2::facet_grid(paste(if(rowVariable != "None") switch(rowVariable, validation = "validation", datasetName = "dataset", classificationName = "analysis"), "~", if(columnVariable != "None") switch(columnVariable, validation = "validation", datasetName = "dataset", classificationName = "analysis")))
  
  if(plot == TRUE)
    print(selectionPlot)
  
  selectionPlot
})