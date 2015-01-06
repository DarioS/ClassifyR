\name{selectionPlot}
\alias{selectionPlot}
\alias{selectionPlot,list-method}
\title{Plot Pair-wise Overlap of Selected Features}
\description{The average pair-wise overlap is computed for every pair of cross-validations. The overlap is converted to a percentage and plotted as a set of boxplots.
       }
\usage{
  \S4method{selectionPlot}{list}(results,
                   comparison = c("within", "classificationName", "validation", "datasetName"),
                   xVariable = c("classificationName", "datasetName", "validation"),
                   boxFillColouring = c("classificationName", "datasetName", "validation", "None"),
                   boxFillColours = NULL,
                   boxLineColouring = c("validation", "classificationName", "datasetName", "None"),
                   boxLineColours = NULL,
                   rowVariable = c("None", "validation", "datasetName", "classificationName"),
                   columnVariable = c("datasetName", "classificationName", "validation", "None"),
                   yMax = 100, fontSizes = c(24, 16, 12), title = "Feature Selection Stability",
                   xLabel = "Analysis", yLabel = "Average Pairwise Common Features (\%)", margin = grid::unit(c(0, 1, 1, 0), "lines"), rotate90 = FALSE, plot = TRUE, parallelParams = bpparam())
}
\arguments{
  \item{results}{A list of \code{\link{ClassifyResult}} objects.}
  \item{comparison}{The aspect of the experimental design to compare. See \code{Details} section for a
                    detailed description.}  
  \item{xVariable}{The factor to make separate boxes in the boxplot for.}
  \item{boxFillColouring}{A factor to colour the boxes by.}
  \item{boxFillColours}{A vector of colours, one for each level of \code{boxFillColouring}. If \code{NULL},
                     a default palette is used.}
  \item{boxLineColouring}{A factor to colour the box lines by.}
  \item{boxLineColours}{A vector of colours, one for each level of \code{boxLineColouring}. If \code{NULL},
                     a default palette is used.}  
  \item{rowVariable}{The slot name that different levels of are plotted as separate rows of boxplots.}     
  \item{columnVariable}{The slot name that different levels of are plotted as separate columns of boxplots.}
  \item{yMax}{The maximum value of the percentage to plot.}
  \item{fontSizes}{A vector of length 3. The first number is the size of the title.
                   The second number is the size of the axes titles. The third number is
                   the size of the axes values.}
  \item{title}{An overall title for the plot.}
  \item{xLabel}{Label to be used for the x-axis.}
  \item{yLabel}{Label to be used for the y-axis of overlap percentages.}
  \item{margin}{The margin to have around the plot.}
  \item{rotate90}{Logical. IF TRUE, the boxplot is horizontal.}
  \item{plot}{Logical. IF TRUE, a plot is produced on the current graphics device.}
  \item{parallelParams}{An object of class \code{\link{MulticoreParam}} or \code{\link{SnowParam}}.}
}
\details{
  Possible values for slot names are \code{"datasetName"}, \code{"classificationName"}, and \code{"validation"}.
  If \code{"None"}, then that graphic element is not used.
  
  Calculating all pair-wise set overlaps can be time-consuming. This stage can be done on multiple CPUs by
  providing the relevant options to \code{parallelParams}.
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
  predicted[, "label"] <- sample(predicted[, "label"])
  result2 <- ClassifyResult("Example", "Differential Variability", LETTERS[1:10], LETTERS[10:1], list(1:100, c(1:5, 11:105)), list(c(1:3), c(4:6), c(1, 6, 7, 9), c(5:8), c(1, 5, 10)),
                            list(predicted), actual, validation = list("leave", 1))
  # selectionPlot(list(result1, result2), xVariable = "classificationName", xLabel = "Analysis", columnVariable = "None", rowVariable = "None", boxColouring = "classificationName") # Wait for namespace problems to be fixed.
}