\name{performancePlot}
\alias{performancePlot}
\alias{performancePlot,list-method}
\title{Plot Performance Measures for Various Classifications}
\description{Draws a graphical summary of a particular performance measure for a list of classifications
       }
\usage{
  \S4method{performancePlot}{list}(results,
                   aggregate = character(),
                   xVariable = c("classificationName", "datasetName", "selectionName", "validation"),
                   performanceName = NULL,
                   boxFillColouring = c("classificationName", "datasetName", "selectionName", "validation", "None"),
                   boxFillColours = NULL,
                   boxLineColouring = c("classificationName", "datasetName", "selectionName", "validation", "None"),
                   boxLineColours = NULL,
                   rowVariable = c("None", "validation", "datasetName", "classificationName", "selectionName"),
                   columnVariable = c("datasetName", "classificationName", "validation", "selectionName", "None"),
                   yMax = 1, fontSizes = c(24, 16, 12, 12), title = NULL,
                   xLabel = "Analysis", yLabel = performanceName,
                   margin = grid::unit(c(0, 0, 0, 0), "lines"), rotate90 = FALSE, showLegend = TRUE, plot = TRUE)
}
\arguments{
  \item{results}{A list of \code{\link{ClassifyResult}} objects.}
  \item{aggregate}{A character vector of the levels of \code{xVariable} to aggregate to a single number
                   by taking the mean. This is partciularly meaningful when the cross-validation is leave-k-out,
                   when k is small.}
  \item{xVariable}{The factor to make separate boxes for.}
  \item{performanceName}{The name of the performance measure to make comparisons of. This is
                         one of the names printed in the Performance Measures field when a
                         \code{\link{ClassifyResult}} object is printed.}
  \item{boxFillColouring}{A factor to colour the boxes by.}
  \item{boxFillColours}{A vector of colours, one for each level of \code{boxFillColouring}.}
  \item{boxLineColouring}{A factor to colour the box lines by.}
  \item{boxLineColours}{A vector of colours, one for each level of \code{boxLineColouring}.}
  \item{rowVariable}{The slot name that different levels of are plotted as separate rows of boxplots.}     
  \item{columnVariable}{The slot name that different levels of are plotted as separate columns of boxplots.}
  \item{yMax}{The maximum value of the percentage to plot.}
  \item{fontSizes}{A vector of length 4. The first number is the size of the title.
                   The second number is the size of the axes titles. The third number is
                   the size of the axes values. The fourth number is the font size of the
                   titles of grouped plots, if any are produced. In other words, when
                   \code{rowVariable} or \code{columnVariable} are not \code{NULL}.}
  \item{title}{An overall title for the plot.}
  \item{xLabel}{Label to be used for the x-axis.}
  \item{yLabel}{Label to be used for the y-axis of overlap percentages.}
  \item{margin}{The margin to have around the plot.}
  \item{rotate90}{Logical. IF \code{TRUE}, the plot is horizontal.}
  \item{showLegend}{If \code{TRUE}, a legend is plotted next to the plot. If FALSE, it is hidden.}
  \item{plot}{Logical. IF \code{TRUE}, a plot is produced on the current graphics device.}
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
  predicted <- list(data.frame(sample = sample(10, 20, replace = TRUE),
                              label = rep(c("Healthy", "Cancer"), each = 10)),
                    data.frame(sample = sample(10, 20, replace = TRUE),
                               label = rep(c("Healthy", "Cancer"), each = 10)),
                    data.frame(sample = sample(10, 20, replace = TRUE),
                               label = rep(c("Healthy", "Cancer"), each = 10)),
                    data.frame(sample = sample(10, 20, replace = TRUE),
                               label = rep(c("Healthy", "Cancer"), each = 10)))
  actual <- factor(rep(c("Healthy", "Cancer"), each = 5))
  result1 <- ClassifyResult("Example", "Differential Expression", "t-test", LETTERS[1:10], LETTERS[10:1], list(1:100, c(1:9, 11:101)), list(c(1:3), c(2, 5, 6), c(1:4), c(5:8), 1:5),
                            predicted, actual, list("fold", 2, 2))
  result1 <- calcPerformance(result1, "f")
  predicted <- data.frame(sample = sample(10, 100, replace = TRUE),
                          label = rep(c("Healthy", "Cancer"), each = 50))
  result2 <- ClassifyResult("Example", "Differential Variability", "F-test", LETTERS[1:10], LETTERS[10:1], list(1:100, c(1:5, 11:105)), list(c(1:3), c(4:6), c(1, 6, 7, 9), c(5:8), c(1, 5, 10)),
                            list(predicted), actual, validation = list("leave", 1))
  result2 <- calcPerformance(result2, "f")                            
  performancePlot(list(result1, result2), performanceName = "Precision-Recall F measure", title = "Comparison", boxLineColouring = "None")
}