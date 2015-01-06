setGeneric("rankPlot", function(results, ...)
{standardGeneric("rankPlot")})

setMethod("rankPlot", "list", 
          function(results, topRanked = seq(10, 100, 10),
                   comparison = c("within", "classificationName", "validation", "datasetName"),
                   lineColourVariable = c("validation", "datasetName", "classificationName", "None"),
                   lineColours = NULL, lineWidth = 1,
                   pointTypeVariable = c("datasetName", "classificationName", "validation", "None"),
                   rowVariable = c("None", "datasetName", "classificationName", "validation"),
                   columnVariable = c("classificationName", "datasetName", "validation", "None"),
                   yMax = 100, fontSizes = c(24, 16, 12, 12, 12), title = "Feature Ranking Stability",
                   xLabelPositions = seq(10, 100, 10), yLabel = "Average Pairwise Common Features (%)",
                   plot = TRUE, parallelParams = bpparam())
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")       
  
  comparison <- match.arg(comparison)
  lineColourVariable <- match.arg(lineColourVariable)
  pointTypeVariable <- match.arg(pointTypeVariable)
  rowVariable <- match.arg(rowVariable)
  columnVariable <- match.arg(columnVariable)
  
  results <- lapply(results, function(result) # Unlist the folds of Resample and Fold.
  {
    if(is.list(result@rankedFeatures[[1]]))
      result@rankedFeatures <- unlist(result@rankedFeatures, recursive = FALSE)
    result
  })
  
  if(comparison == "within")
  {
    plotData <- do.call(rbind, bplapply(results, function(result)
    {
      rankedFeatures <- result@rankedFeatures
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
                 validation = rep(validationText, length(topRanked)),
                 top = topRanked,
                 overlap = averageOverlap)
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

    plotData <- do.call(rbind, unname(unlist(bplapply(compareIndices, function(indicesSet)
    {
      if(length(indicesSet) == 2)
        indiciesCombinations <- list(indicesSet)
      else
        indiciesCombinations <- lapply(1:length(indicesSet), function(index) c(indicesSet[index], indicesSet[-index]))
      unname(lapply(indiciesCombinations, function(indiciesCombination)
      {
        aDataset <- results[[indiciesCombination[1]]]
        otherDatasets <- results[indiciesCombination[-1]]   
        averageOverlap <- rowMeans(do.call(cbind, lapply(aDataset@rankedFeatures, function(features)
        {
          do.call(cbind, unname(lapply(otherDatasets, function(anotherDataset)
          {
            sapply(anotherDataset@rankedFeatures, function(otherFeatures)
            {
              sapply(topRanked, function(top)
              {
                length(intersect(features[1:top], otherFeatures[1:top])) / top * 100
              })          
            })
          })))
        })))
        validationText <- .validationText(aDataset)
        
        data.frame(dataset = rep(aDataset@datasetName, length(topRanked)),
                   analysis = rep(aDataset@classificationName, length(topRanked)),
                   validation = rep(validationText, length(topRanked)),
                   top = topRanked,
                   overlap = averageOverlap)
      }))
    }, BPPARAM = parallelParams), recursive = FALSE)))
  }
  
  plotData[, "dataset"] <- factor(plotData[, "dataset"], levels = unique(sapply(results, function(result) result@datasetName)))
  plotData[, "analysis"] <- factor(plotData[, "analysis"], levels = unique(sapply(results, function(result) result@classificationName)))
  plotData[, "validation"] <- factor(plotData[, "validation"], levels = unique(plotData[, "validation"])) # Order in which validations were provided.
  
  overlapPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = top, y = overlap,
                          colour = switch(lineColourVariable, validation = validation, datasetName = dataset, classificationName = analysis),
                          shape = switch(pointTypeVariable, validation = validation, datasetName = dataset, classificationName = analysis)), environment = environment()) +
                          ggplot2::geom_line(size = lineWidth) + ggplot2::geom_point(size = 3) + ggplot2::scale_x_continuous(breaks = xLabelPositions) + ggplot2::scale_y_continuous(limits = c(0, yMax)) + ggplot2::xlab("Top Features") + ggplot2::ylab(yLabel) +
                          ggplot2::ggtitle(title) + ggplot2::labs(colour = switch(lineColourVariable, validation = "Validation", datasetName = "Dataset", classificationName = "Analysis"), shape = switch(pointTypeVariable, validation = "Validation", datasetName = "Dataset", classificationName = "Analysis")) +
                          ggplot2::theme(axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), legend.title = ggplot2::element_text(size = fontSizes[4]), legend.text = ggplot2::element_text(size = fontSizes[5]), plot.title = ggplot2::element_text(size = fontSizes[1]), plot.margin = grid::unit(c(0, 0, 1, 0), "lines"), legend.margin = grid::unit(-1, "lines"))
  
  if(rowVariable != "None" || columnVariable != "None")
    overlapPlot <- overlapPlot + ggplot2::facet_grid(paste(if(rowVariable != "None") switch(rowVariable, validation = "validation", datasetName = "dataset", classificationName = "analysis"), "~", if(columnVariable != "None") switch(columnVariable, validation = "validation", datasetName = "dataset", classificationName = "analysis")))
  
  if(plot == TRUE)
    print(overlapPlot)
  
  overlapPlot
})
