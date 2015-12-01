\name{ROCplot}
\alias{ROCplot}
\alias{ROCplot,list-method}
\title{Plot Receiver Operating Curve Graphs for Classification Results}
\description{The average pair-wise overlap is computed for every pair of cross-validations. The overlap is converted to a percentage and plotted as lineplots.
       }
\usage{
  \S4method{ROCplot}{list}(results, nBins = sapply(results, totalPredictions),
                   lineColourVariable = c("classificationName", "datasetName", "selectionName", "validation", "None"), lineColours = NULL,
                   lineWidth = 1, fontSizes = c(24, 16, 12, 12, 12), labelPositions = seq(0.0, 1.0, 0.2),
                   plotTitle = "ROC", legendTitle = NULL, xLabel = "False Positive Rate", yLabel = "True Positive Rate",
                   plot = TRUE, showAUC = TRUE)
}
\arguments{
  \item{results}{A list of \code{\link{ClassifyResult}} objects.}
  \item{nBins}{The number of intervals to group the samples' scores into. By default, there are as
               many bins as there were predictions made, for each result object.}
  \item{lineColourVariable}{The slot name that different levels of are plotted as
                            different line colours.}
  \item{lineColours}{A vector of colours for different levels of the line colouring parameter. If \code{NULL},
                     a default palette is used.}       
  \item{lineWidth}{A single number controlling the thickness of lines drawn.}
  \item{fontSizes}{A vector of length 5. The first number is the size of the title.
                   The second number is the size of the axes titles and AUC text, if it
                   is not part of the legend. The third number is the size of the axes values. The
                   fourth number is the size of the legends' titles. The fifth number is
                   the font size of the legend labels.}
  \item{labelPositions}{Locations where to put labels on the x and y axes.}
  \item{plotTitle}{An overall title for the plot.}
  \item{legendTitle}{A default name is used if the value is \code{NULL}. Otherwise a character name can
                   be provided.}
  \item{xLabel}{Label to be used for the x-axis of false positive rate.}
  \item{yLabel}{Label to be used for the y-axis of true positive rate.}
  \item{plot}{Logical. If \code{TRUE}, a plot is produced on the current graphics device.}
  \item{showAUC}{Logical. If \code{TRUE}, the AUC value of each result is added to its legend text.}
}
\details{
  Possible values for slot names are \code{"datasetName"}, \code{"classificationName"}, and
  \code{"validation"}. If \code{"None"}, then any lines drawn will be black.
  
  The scores stored in the results should be higher if the sample is more likely to be from the second class,
  based on the levels of the actual classes. The scores must be in a column named "score".
  
  For cross-validated classification, all predictions from all iterations are considered simultaneously,
  to calculate one curve per classification.
  
  The number of bins determines how many pairs of TPR and FPR points will be used to draw the plot.
  A higher number will result in a smoother ROC curve.
  
  The AUC is calculated using the trapezoidal rule.
}
\value{
  An object of class \code{ggplot} and a plot on the current graphics device, if \code{plot} is \code{TRUE}.
}
\author{Dario Strbenac}

\examples{
  predicted <- list(data.frame(sample = c(1, 8, 15, 3, 11, 20, 19, 18), score = c(0.11, 0.32, 0.47, 0.24, 0.87, 0.80, 0.40, 0.75)),
                    data.frame(sample = c(11, 18, 15, 4, 6, 10, 11, 12), score = c(0.55, 0.44, 0.67, 0.44, 0.67, 0.80, 0.40, 0.60)))
  actual <- factor(c(rep("Healthy", 10), rep("Cancer", 10)), levels = c("Healthy", "Cancer"))
  result1 <- ClassifyResult("Example", "Differential Expression", "t-test", LETTERS[1:10], LETTERS[10:1], list(1:100, c(1:9, 11:101)), list(sample(10, 10), sample(10, 10)),
                            predicted, actual, list("fold", 2, 1))
  predicted[[1]][, "score"][c(2, 6)] <- c(0.60, 0.40)
  result2 <- ClassifyResult("Example", "Differential Variability", "F-test", LETTERS[1:10], LETTERS[10:1], list(1:100, c(1:5, 11:105)), list(sample(10, 10), sample(10, 10)),
                            predicted, actual, validation = list("fold", 2, 1))
  ROCplot(list(result1, result2), lineColourVariable = "classificationName", plotTitle = "Ovarian Cancer ROC")
}