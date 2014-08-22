setGeneric("plotFeatureClasses", function(expression, ...)
           {standardGeneric("plotFeatureClasses")})

setMethod("plotFeatureClasses", "matrix", 
          function(expression, classes, ...)
{
  groupsTable <- data.frame(class = classes)
  features <- rownames(expression)
  rownames(expression) <- NULL
  rownames(groupsTable) <- colnames(expression)
  exprSet <- ExpressionSet(expression, AnnotatedDataFrame(groupsTable))
  if(length(features) > 0) featureNames(exprSet) <- features
  plotFeatureClasses(exprSet, ...)
})

setMethod("plotFeatureClasses", "ExpressionSet", 
          function(expression, rows, assayType)
{
  classes <- pData(expression)[, "class"]
  features <- featureNames(expression)
  expression <- exprs(expression)
  invisible(lapply(rows, function(featureRow)
  {
    plotData <- data.frame(expr = expression[featureRow, ], classes)
    densPlot <- ggplot(plotData, aes(x = expr, colour = classes)) +
                stat_density(aes(y = ..density..), geom = "path", position = "identity", size = 1) +
                scale_colour_manual("Vital", values = c("red", "blue")) +
                theme(axis.title.x = element_blank(), axis.line = element_blank()) 
    if(assayType == "array")
      densPlot <- densPlot + scale_x_continuous(limits = c(2, 16))
    stripPlot <- ggplot(plotData, aes(x = classes, y = expr)) +
                 geom_dotplot(dotsize = 0.5, binaxis = "y", stackdir = "center", position = "dodge", aes(colour = classes)) +
                 scale_colour_manual("Vital", values = c("red", "blue"))
      if(assayType == "array")
        stripPlot <- stripPlot + ylab(expression(log[2](expression)))
    else
      stripPlot <- stripPlot + ylab("Transcripts Per Million (T. P. M.)")
    if(assayType == "array")
      stripPlot <- stripPlot + scale_y_continuous(limits = c(2, 16))
    stripPlot <- stripPlot + coord_flip()
    if(is.numeric(rows)) featureName <- features[featureRow] else featureName <-featureRow
    grid.arrange(densPlot, stripPlot, nrow = 2,
                 main = textGrob(featureName, gp = gpar(fontsize = 18), vjust = 1))
  }))
})