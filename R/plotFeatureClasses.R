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
          function(expression, rows, plot = c("both", "density", "stripchart"),
                   expressionLabel = expression(log[2](expression)), expressionLimits = c(2, 16),
                   colours = c("blue", "red"))
{
  plot <- match.arg(plot)
  classes <- pData(expression)[, "class"]
  features <- featureNames(expression)
  expression <- exprs(expression)
  invisible(lapply(rows, function(featureRow)
  {
    plotData <- data.frame(expr = expression[featureRow, ], classes)
    if(plot %in% c("both", "density"))
    {
      densPlot <- ggplot(plotData, aes(x = expr, colour = classes)) +
                  stat_density(aes(y = ..density..), geom = "path", position = "identity", size = 1) +
                  scale_colour_manual("Vital", values = c("red", "blue")) +
                  theme(axis.title.x = element_blank(), axis.line = element_blank()) +
                  scale_x_continuous(limits = expressionLimits)
    }
    if(plot %in% c("both", "stripchart"))
    {
      stripPlot <- ggplot(plotData, aes(x = classes, y = expr)) +
                   geom_dotplot(dotsize = 0.5, binaxis = "y", stackdir = "center", position = "dodge", aes(colour = classes)) +
                   scale_colour_manual("Vital", values = c("red", "blue")) + ylab(expressionLabel) + scale_y_continuous(limits = expressionLimits)
      stripPlot <- stripPlot + coord_flip()
    }
    if(is.numeric(rows)) featureName <- features[featureRow] else featureName <-featureRow
    if(plot == "both")
       grid.arrange(densPlot, stripPlot, nrow = 2,
                    main = textGrob(featureName, gp = gpar(fontsize = 18), vjust = 1))
    else if(plot == "density")
      densPlot + ggtitle(featureName)
    else
      stripPlot + ggtitle(featureName)
  }))
})