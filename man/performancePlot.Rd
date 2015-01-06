\name{performancePlot}
\alias{performancePlot}
\alias{performancePlot,list-method}
\title{Plot Performance Measures for Various Classifications}
\description{Draws a graphical summary of a particular performance measure for a list of classifications
       }
\usage{
  \S4method{performancePlot}{list}(results,
                   aggregate = FALSE,
                   xVariable = c("classificationName", "datasetName", "validation"),
                   performanceName = NULL,
                   boxFillColouring = c("classificationName", "datasetName", "validation", "None"),
                   boxFillColours = NULL,
                   boxLineColouring = c("classificationName", "datasetName", "validation", "None"),
                   boxLineColours = NULL,
                   yMax = 1, fontSizes = c(24, 16, 12), title = NULL,
                   xLabel = "Classification", yLabel = performanceName,
                   margin = grid::unit(c(0, 1, 1, 0), "lines"), rotate90 = FALSE, plot = TRUE)
}
\arguments{
  \item{results}{A list of \code{\link{ClassifyResult}} objects.}
  \item{aggregate}{If TRUE, and there are multiple values of a particular performance measure stored,
                   such as when bootstrap resampling is done, the performance measure is summarised
                   to a single number by considering all predictions simultaneously.}
  \item{xVariable}{The factor to make separate boxes for.}
  \item{performanceName}{The name of the performance measure to make comparisons of. This is
                         one of the names printed in the Performance Measures field when a
                         \code{\link{ClassifyResult}} object is printed.}
  \item{boxFillColouring}{A factor to colour the boxes by.}
  \item{boxFillColours}{A vector of colours, one for each level of \code{boxFillColouring}.}
  \item{boxLineColouring}{A factor to colour the box lines by.}
  \item{boxLineColours}{A vector of colours, one for each level of \code{boxLineColouring}.}  
  \item{yMax}{The maximum value of the percentage to plot.}
  \item{fontSizes}{A vector of length 3. The first number is the size of the title.
                   The second number is the size of the axes titles. The third number is
                   the size of the axes values.}
  \item{title}{An overall title for the plot.}
  \item{xLabel}{Label to be used for the x-axis.}
  \item{yLabel}{Label to be used for the y-axis of overlap percentages.}
  \item{margin}{The margin to have around the plot.}
  \item{rotate90}{Logical. IF TRUE, the plot is horizontal.}
  \item{plot}{Logical. IF TRUE, a plot is produced on the current graphics device.}
}
\details{
  Possible values for slot names are \code{"datasetName"}, \code{"classificationName"}, and \code{"validation"}.
  If \code{"None"}, then that graphic element is not used.
  
  If there are multiple values for a performance measure in a single result object, it is plotted as a boxplot,
  unless \code{aggregate} is \code{TRUE}, in which case the all predictions in a single result object are considered simultaneously, so that only one performance number is calculated, and a barchart is plotted.
}
\value{
  An object of class \code{ggplot} and a plot on the current graphics device, if \code{plot} is \code{TRUE}.
}
\author{Dario Strbenac}

\examples{
  predicted <- data.frame(sample = sample(10, 100, replace = TRUE),
                          label = rep(c("Healthy", "Cancer"), each = 50))
  actual <- factor(rep(c("Healthy", "Cancer"), each = 5))
  result1 <- ClassifyResult("Example", "Differential Expression", LETTERS[1:10], LETTERS[10:1], list(1:100, c(1:9, 11:101)), list(c(1:3), c(2, 5, 6), c(1:4), c(5:8), 1:5),
                            list(predicted), actual, list("fold", 100, 5))
  result1 <- calcPerformance(result1, "f")
  predicted[, "label"] <- sample(predicted[, "label"])
  result2 <- ClassifyResult("Example", "Differential Variability", LETTERS[1:10], LETTERS[10:1], list(1:100, c(1:5, 11:105)), list(c(1:3), c(4:6), c(1, 6, 7, 9), c(5:8), c(1, 5, 10)),
                            list(predicted), actual, validation = list("leave", 1))
  result2 <- calcPerformance(result2, "f")                            
  performancePlot(list(result1, result2), performanceName = "Precision-Recall F measure", title = "Comparison")
}