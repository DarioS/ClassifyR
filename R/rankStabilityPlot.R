setGeneric("rankStabilityPlot", function(results, ...)
{standardGeneric("rankStabilityPlot")})

setMethod("rankStabilityPlot", "list", 
          function(results, topRanked = seq(10, 100, 10),
                   lineColourVariable = c("validation", "datasetName", "classificationName", "None"),
                   pointTypeVariable = c("datasetName", "classificationName", "validation", "None"),
                   rowVariable = c("None", "datasetName", "classificationName", "validation"),
                   columnVariable = c("classificationName", "datasetName", "validation", "None"),
                   fontSizes = c(24, 16, 12, 12, 12), plot = TRUE, parallelParams = bpparam())
{
  lineColourVariable <- match.arg(lineColourVariable)
  pointTypeVariable <- match.arg(pointTypeVariable)
  rowVariable <- match.arg(rowVariable)
  columnVariable <- match.arg(columnVariable)
            
  plotData <- do.call(rbind, lapply(results, function(result)
  {
    if(is.list(result@rankedFeatures[[1]]))
      rankedFeatures <- unlist(result@rankedFeatures, recursive = FALSE)
    else
      rankedFeatures <- result@rankedFeatures
    averageOverlap <- rowMeans(do.call(cbind, bpmapply(function(features, index)
    {
      otherFeatures <- rankedFeatures[(index + 1):length(rankedFeatures)]
      sapply(otherFeatures, function(other)
      {
        sapply(topRanked, function(top)
        {
          
          length(intersect(features[1:top], other[1:top])) / top * 100
        })
      })
    }, rankedFeatures[1:(length(rankedFeatures) - 1)], 1:(length(rankedFeatures) - 1), SIMPLIFY = FALSE, BPPARAM = parallelParams)))
    
    validationText <- if(result@validation[[1]] == "fold") "Resample and Fold"
                      else if(result@validation[[1]] == "split") "Resample and Split"
                      else "Leave Out"
    data.frame(dataset = rep(result@datasetName, length(topRanked)),
               analysis = rep(result@classificationName, length(topRanked)),
               validation = rep(validationText, length(topRanked)),
               top = topRanked,
               overlap = averageOverlap)
  }))

  stabilityPlot <- ggplot2::ggplot(plotData, ggplot2::aes(x = top, y = overlap,
                          colour = switch(lineColourVariable, validation = validation, datasetName = dataset, classificationName = analysis),
                          shape = switch(pointTypeVariable, validation = validation, datasetName = dataset, classificationName = analysis)), environment = environment()) +
                          ggplot2::geom_line() + ggplot2::geom_point(size = 3) + ggplot2::scale_x_continuous(breaks = topRanked) + ggplot2::scale_y_continuous(limits = c(0, 100)) + ggplot2::xlab("Top Features") + ggplot2::ylab("Average Pairwise Common Features (%)") +
                          ggplot2::ggtitle("Feature Selection Stability") + ggplot2::labs(colour = switch(lineColourVariable, validation = "Validation", datasetName = "Dataset", classificationName = "Classification"), shape = switch(pointTypeVariable, validation = "Validation", datasetName = "Dataset", classificationName = "Classification")) +
                          ggplot2::theme(axis.title = ggplot2::element_text(size = fontSizes[2]), axis.text = ggplot2::element_text(colour = "black", size = fontSizes[3]), legend.title = ggplot2::element_text(size = fontSizes[4]), legend.text = ggplot2::element_text(size = fontSizes[5]), plot.title = ggplot2::element_text(size = fontSizes[1]))
  
  if(rowVariable != "None" || columnVariable != "None")
    stabilityPlot <- stabilityPlot + ggplot2::facet_grid(paste(if(rowVariable != "None") switch(rowVariable, validation = "validation", datasetName = "dataset", classificationName = "analysis"), "~", if(columnVariable != "None") switch(columnVariable, validation = "validation", datasetName = "dataset", classificationName = "analysis")))
  
  if(plot == TRUE)
    print(stabilityPlot)
  
  stabilityPlot
})
