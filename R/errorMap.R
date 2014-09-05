setGeneric("errorMap", function(results, ...)
           {standardGeneric("errorMap")})

setMethod("errorMap", "list", 
          function(results,
                   errorColours = c("#FFFFFF", "#E0E0E0", "#BABABA", "#888888", "#000000"),
                   classColours = c("blue", "red"), fontSizes = c(24, 16, 12, 12, 12))
{
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
  
  classesPlot <- ggplot(data.frame(Class = knownClasses), aes(Class, factor(1))) +
                 scale_fill_manual(values = classColours) + geom_tile(aes(fill = Class, height = 10)) +
                 scale_x_discrete(expand = c(0, 0), breaks = NULL) +
                 scale_y_discrete(expand = c(0, 0), breaks = NULL) +
                 labs(x = '', y = '') + theme(plot.margin = unit(c(2, 1, 0, 1), "cm"),
                                              legend.title = element_text(size = fontSizes[4]),
                                              legend.text = element_text(size = fontSizes[5]))
  
  errorPlot <- ggplot(plotData, aes(name, type)) + geom_tile(aes(fill = Error)) +
               scale_fill_manual(values = errorColours, drop = FALSE) + scale_x_discrete(expand = c(0, 0)) +
               scale_y_discrete(expand = c(0, 0)) + theme_bw() +
               theme(axis.ticks = element_blank(),
               axis.text.x = element_text(angle = 45, hjust = 1, size = fontSizes[3]),
               axis.text.y = element_text(size = fontSizes[3]),
               axis.title = element_text(size = fontSizes[2]),
               plot.margin = unit(c(0, 1, 1, 1), "cm"),
               legend.title = element_text(size = fontSizes[4]),
               legend.text = element_text(size = fontSizes[5])) + labs(x = "Sample Name", y = "Result")
  
  classGrob <- ggplot_gtable(ggplot_build(classesPlot))
  errorGrob <- ggplot_gtable(ggplot_build(errorPlot))
  commonWidth <- unit.pmax(classGrob[["widths"]], errorGrob[["widths"]])
  classGrob[["widths"]] <- commonWidth
  errorGrob[["widths"]] <- commonWidth
  
  grid.arrange(classGrob, errorGrob, nrow = 2, heights = c(length(results), 4 * length(results)),
               main = textGrob("Error Comparison", gp = gpar(fontsize = fontSizes[1]), vjust = 1))
})