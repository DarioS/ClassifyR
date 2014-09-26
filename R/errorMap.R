setGeneric("errorMap", function(results, ...)
           {standardGeneric("errorMap")})

setMethod("errorMap", "list", 
          function(results,
                   errorColours = c("#FFFFFF", "#E0E0E0", "#BABABA", "#888888", "#000000"),
                   classColours = c("blue", "red"), fontSizes = c(24, 16, 12, 12, 12))
{
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("The package 'ggplot2' could not be found. Please install it.")  
  if(!requireNamespace("grid", quietly = TRUE))
    stop("The package 'grid' could not be found. Please install it.")
  if(!requireNamespace("gridExtra", quietly = TRUE))
    stop("The package 'gridExtra' could not be found. Please install it.")       
            
  errorBinEnds <- seq(0, 1, 1/length(errorColours))
  errors <- lapply(results, function(result)
  {
    resultTable <- do.call(rbind, predictions(result))
    sampleErrors <- by(resultTable, resultTable[, "sample"],
                      function(samplePredictions)
                        sum(samplePredictions[, "predicted"] != result@actualClasses[samplePredictions[1, "sample"]]))
    cut(sampleErrors / table(resultTable[, "sample"]), errorBinEnds, include.lowest = TRUE)
  })
  
  knownClasses <- actualClasses(results[[1]])
  
  ordering <- order(knownClasses)
  knownClasses <- knownClasses[ordering]
  errors <- lapply(errors, function(resultErrors) resultErrors[ordering])
  
  plotData <- data.frame(name = rep(sampleNames(results[[1]])[ordering], length(results)),
                         type = rep(names(results), sapply(errors, length)),
                         class = rep(knownClasses, length(results)),
                         Error = unlist(errors))
  
  classesPlot <- ggplot2::ggplot(data.frame(Class = knownClasses), ggplot2::aes(1:length(knownClasses), factor(1)), environment = environment()) +
                 ggplot2::scale_fill_manual(values = classColours) + ggplot2::geom_tile(aes(fill = Class, height = 10)) +
                 ggplot2::scale_x_discrete(expand = c(0, 0), breaks = NULL, limits = c(1, length(knownClasses))) +
                 ggplot2::scale_y_discrete(expand = c(0, 0), breaks = NULL) +
                 ggplot2::labs(x = '', y = '') + ggplot2::theme(plot.margin = grid::unit(c(2, 1, 0, 1), "cm"),
                                                                legend.title = ggplot2::element_text(size = fontSizes[4]),
                                                                legend.text = ggplot2::element_text(size = fontSizes[5]))
  
  errorPlot <- ggplot2::ggplot(plotData, aes(name, type)) + ggplot2::geom_tile(ggplot2::aes(fill = Error)) +
               ggplot2::scale_fill_manual(values = errorColours, drop = FALSE) + ggplot2::scale_x_discrete(expand = c(0, 0)) +
               ggplot2::scale_y_discrete(expand = c(0, 0)) + ggplot2::theme_bw() +
               ggplot2::theme(axis.ticks = ggplot2::element_blank(),
               axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = fontSizes[3]),
               axis.text.y = ggplot2::element_text(size = fontSizes[3]),
               axis.title = ggplot2::element_text(size = fontSizes[2]),
               plot.margin = grid::unit(c(0, 1, 1, 1), "cm"),
               legend.title = ggplot2::element_text(size = fontSizes[4]),
               legend.text = ggplot2::element_text(size = fontSizes[5])) + ggplot2::labs(x = "Sample Name", y = "Result")
  
  classGrob <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(classesPlot))
  errorGrob <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(errorPlot))
  commonWidth <- grid::unit.pmax(classGrob[["widths"]], errorGrob[["widths"]])
  classGrob[["widths"]] <- commonWidth
  errorGrob[["widths"]] <- commonWidth
  
  gridExtra::grid.arrange(classGrob, errorGrob, nrow = 2, heights = c(length(results), 4 * length(results)),
               main = textGrob("Error Comparison", gp = gpar(fontsize = fontSizes[1]), vjust = 1))
})