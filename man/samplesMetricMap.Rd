\name{samplesMetricMap}
\alias{samplesMetricMap}
\alias{samplesMetricMap,list-method}
\alias{samplesMetricMap,matrix-method}
\title{Plot a Grid of Sample Error Rates or Accuracies}
\description{A grid of coloured tiles is drawn. There is one column for each sample
and one row for each classification result.
            }
\usage{
  \S4method{samplesMetricMap}{list}(results,
                     comparison = c("classificationName", "datasetName", "selectionName",
                                "validation"),
                 metric = c("error", "accuracy"), featureValues = NULL, featureName = NULL,
                     metricColours = list(c("#3F48CC", "#6F75D8", "#9FA3E5", "#CFD1F2", "#FFFFFF"),
                                        c("#880015", "#A53F4F", "#C37F8A", "#E1BFC4", "#FFFFFF")),
                     classColours = c("#3F48CC", "#880015"), groupColours = c("darkgreen", "yellow2"),
                     fontSizes = c(24, 16, 12, 12, 12),
                  mapHeight = 4, title = "Error Comparison", showLegends = TRUE,
                 xAxisLabel = "Sample Name", showXtickLabels = TRUE,
                 yAxisLabel = "Analysis", showYtickLabels = TRUE,
                 legendSize = grid::unit(1, "lines"), plot = TRUE)
  \S4method{samplesMetricMap}{matrix}(results, classes,
                     metric = c("error", "accuracy"),
                     featureValues = NULL, featureName = NULL,
                     metricColours = list(c("#3F48CC", "#6F75D8", "#9FA3E5", "#CFD1F2", "#FFFFFF"),
                                        c("#880015", "#A53F4F", "#C37F8A", "#E1BFC4", "#FFFFFF")),
                     classColours = c("#3F48CC", "#880015"), groupColours = c("darkgreen", "yellow2"),
                     fontSizes = c(24, 16, 12, 12, 12),
                  mapHeight = 4, title = "Error Comparison", showLegends = TRUE,
                 xAxisLabel = "Sample Name", showXtickLabels = TRUE,
                 yAxisLabel = "Analysis", showYtickLabels = TRUE,
                 legendSize = grid::unit(1, "lines"), plot = TRUE)                 
}
\arguments{
  \item{results}{A list of \code{\link{ClassifyResult}} objects. Could also be a matrix, for
                 backwards compatibility.}
  \item{classes}{If \code{results} is a matrix, this is a factor vector of the same length
                 as the number of columns that \code{results} has.}                 
  \item{comparison}{The aspect of the experimental design to compare.}
  \item{metric}{The sample-wise metric to plot.}
  \item{featureValues}{If not NULL, can be a named factor or named numeric vector specifying
                some variable of interest to plot underneath the class bar.}
  \item{featureName}{A label describing the information in \code{featureValues}. It must be
                     specified if \code{featureValues} is.}                
  \item{metricColours}{A vector of colours for metric levels.}
  \item{classColours}{Either a vector of colours for class levels if both classes should
                      have same colour, or a list of length 2, with each component
                      being a vector of the same length. The vector has the colour gradient
                      for each class.}
  \item{groupColours}{A vector of colours for group levels. Only useful if \code{groups} is not NULL.}
  \item{fontSizes}{A vector of length 5. The first number is the size of the title.
                   The second number is the size of the axes titles. The third number is
                   the size of the axes values. The fourth number is the size of the
                   legends' titles. The fifth number is the font size of the legend labels.}
  \item{mapHeight}{Height of the map, relative to the height of the class colour bar.}     
  \item{title}{The title to place above the plot.}     
  \item{showLegends}{Logical. IF FALSE, the legend is not drawn.}
  \item{xAxisLabel}{The name plotted for the x-axis. NULL suppresses label.} 
  \item{showXtickLabels}{Logical. IF FALSE, the x-axis labels are hidden.}
  \item{showYtickLabels}{Logical. IF FALSE, the y-axis labels are hidden.}
  \item{yAxisLabel}{The name plotted for the y-axis. NULL suppresses label.}
  \item{legendSize}{The size of the boxes in the legends.}
  \item{plot}{Logical. IF \code{TRUE}, a plot is produced on the current graphics device.}
}
\details{
The names of \code{results} determine the row names that will be in the plot.
The length of \code{metricColours} determines how many bins the metric values will
be discretised to.
}
\value{
  A plot is produced and a grob is returned that can be saved to a graphics device.
}
\author{Dario Strbenac}

\examples{
  predicted <- data.frame(sample = LETTERS[sample(10, 100, replace = TRUE)],
                          class = rep(c("Healthy", "Cancer"), each = 50))
  actual <- factor(rep(c("Healthy", "Cancer"), each = 5), levels = c("Healthy", "Cancer"))
  features <- sapply(1:100, function(index) paste(sample(LETTERS, 3), collapse = ''))
  result1 <- ClassifyResult("Example", "Differential Expression", "t-test",
                            LETTERS[1:10], features, 100, list(1:100), list(sample(10, 10)),
                            list(function(oracle){}), list(predicted), actual,
                            list("permuteFold", 100, 5))
  predicted[, "class"] <- sample(predicted[, "class"])
  result2 <- ClassifyResult("Example", "Differential Variability", "Bartlett Test",
                            LETTERS[1:10], features, 100, list(1:100), list(sample(10, 10)),
                            list(function(oracle){}), list(predicted), actual,
                            validation = list("leave", 2))
  result1 <- calcCVperformance(result1, "sample error")
  result2 <- calcCVperformance(result2, "sample error")
  groups <- factor(rep(c("Male", "Female"), length.out = 10))
  names(groups) <- LETTERS[1:10]
  cholesterol <- c(4.0, 5.5, 3.9, 4.9, 5.7, 7.1, 7.9, 8.0, 8.5, 7.2)
  names(cholesterol) <- LETTERS[1:10]
  
  wholePlot <- samplesMetricMap(list(Gene = result1, Protein = result2))
  wholePlot <- samplesMetricMap(list(Gene = result1, Protein = result2),
                                featureValues = groups, featureName = "Gender")
  wholePlot <- samplesMetricMap(list(Gene = result1, Protein = result2),
                                featureValues = cholesterol, featureName = "Cholesterol")                                
}
