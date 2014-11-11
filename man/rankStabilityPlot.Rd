\name{rankStabilityPlot}
\alias{rankStabilityPlot}
\alias{rankStabilityPlot,list-method}
\title{Plot Pair-wise Overlap of Ranked Features}
\description{The average pair-wise overlap is computed for every pair of cross-validations within a single analysis. The overlap is converted to a percentage and plotted.
       }
\usage{
  \S4method{rankStabilityPlot}{list}(results, topRanked = seq(10, 100, 10),
                              lineColourVariable = c("validation", "datasetName", "classificationName", "None"),
                              pointTypeVariable = c("datasetName", "classificationName", "validation", "None"),
                              rowVariable = c("None", "datasetName", "classificationName", "validation"),
                              columnVariable = c("classificationName", "datasetName", "validation", "None"),
                              fontSizes = c(24, 16, 12, 12, 12), plot = TRUE, parallelParams = bpparam())
}
\arguments{
  \item{results}{A list of \code{\link{ClassifyResult}} objects.}
  \item{topRanked}{A sequence of thresholds of number of the best features to use for overlapping.}
  \item{lineColourVariable}{The slot name that different levels of are plotted as
                            different line colours.}
  \item{pointTypeVariable}{The slot name that different levels of are plotted as
                            different point shapes on the lines.}
  \item{rowVariable}{The slot name that different levels of are plotted as separate rows of plots.}     
  \item{columnVariable}{The slot name that different levels of are plotted as separate columns of plots.}
  \item{fontSizes}{A vector of length 5. The first number is the size of the title.
                   The second number is the size of the axes titles. The third number is
                   the size of the axes values. The fourth number is the size of the
                   legends' titles. The fifth number is the font size of the legend labels.}
  \item{plot}{Logical. IF TRUE, a plot is produced on the current graphics device.}
  \item{parallelParams}{An object of class \code{\link{MulticoreParam}} or \code{\link{SnowParam}}.}
}
\details{
  Possible values for slot names are \code{"datasetName"}, \code{"classificationName"}, and \code{"validation"}. If \code{"None"}, then that graphic element is not used.
  
  Calculating all pair-wise set overlaps can be time-consuming. This stage can be done on multiple CPUs by providing the
  relevant options to \code{parallelParams}.
}
\value{
  An object of class \code{ggplot} and a plot on the current graphics device, if \code{plot} is \code{TRUE}.
}
\author{Dario Strbenac}

\examples{
  predicted <- data.frame(sample = sample(10, 100, replace = TRUE),
                            predicted = rep(c("Healthy", "Cancer"), each = 50))
  actual <- factor(rep(c("Healthy", "Cancer"), each = 5))
  result1 <- ClassifyResult("Example", "Differential Expression", LETTERS[1:10], LETTERS[10:1], list(1:100, c(1:9, 11:101)), list(sample(10, 10), sample(10, 10)),
                            list(predicted), actual, list("fold", 100, 5))
  predicted[, "predicted"] <- sample(predicted[, "predicted"])
  result2 <- ClassifyResult("Example", "Differential Deviation", LETTERS[1:10], LETTERS[10:1], list(1:100, c(1:5, 11:105)), list(sample(10, 10), sample(10, 10)),
                            list(predicted), actual, validation = list("leave", 1))
  # rankStabilityPlot(list(result1, result2), pointTypeVariable = "classificationName") # Wait for namespace problems to be fixed.
}