\name{ROCplot}
\alias{ROCplot}
\alias{ROCplot,list-method}
\title{Plot Receiver Operating Curve Graphs for Classification Results}
\description{Creates one ROC plot or multiple ROC plots for a list of ClassifyR objects. One plot
             is created if the data set has two classes and multiple plots are created if the data
             set has three or more classes.
       }
\usage{
  \S4method{ROCplot}{list}(results, nBins = sapply(results, totalPredictions),
          comparisonVariable = c("classificationName", "datasetName", "selectionName",
        "validation", "None"), lineColours = NULL,
          lineWidth = 1, fontSizes = c(24, 16, 12, 12, 12), labelPositions = seq(0.0, 1.0, 0.2),
        plotTitle = "ROC", legendTitle = NULL, xLabel = "False Positive Rate",
        yLabel = "True Positive Rate", plot = TRUE, showAUC = TRUE)
}
\arguments{
  \item{results}{A list of \code{\link{ClassifyResult}} objects.}
  \item{nBins}{The number of intervals to group the samples' scores into. By default, there are as
               many bins as there were predictions made, for each result object.}
  \item{comparisonVariable}{If the data set has two classes, then the slot name with factor levels
                           to be used for colouring the lines. Possible values are \code{"datasetName"},
                           \code{"classificationName"}, and \code{"validation"}. If \code{"None"},
                           then all lines drawn will be black. Otherwise, it specifies the variable used
                           for plot facetting.}
  \item{lineColours}{A vector of colours for different levels of the line colouring parameter.
                     If \code{NULL}, a default colour palette is used.}
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
  The scores stored in the results should be higher if the sample is more likely to be from the class
  which the score is associated with. The score for each class must be in a column which has a column
  name equal to the class name.
  
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
  predicted <- list(data.frame(sample = LETTERS[c(1, 8, 15, 3, 11, 20, 19, 18)],
                               Healthy = c(0.89, 0.68, 0.53, 0.76, 0.13, 0.20, 0.60, 0.25),
                               Cancer = c(0.11, 0.32, 0.47, 0.24, 0.87, 0.80, 0.40, 0.75)),
                    data.frame(sample = LETTERS[c(11, 18, 15, 4, 6, 10, 11, 12)],
                               Healthy = c(0.45, 0.56, 0.33, 0.56, 0.33, 0.20, 0.60, 0.40),
                               Cancer = c(0.55, 0.44, 0.67, 0.44, 0.67, 0.80, 0.40, 0.60)))
  actual <- factor(c(rep("Healthy", 10), rep("Cancer", 10)), levels = c("Healthy", "Cancer"))
  result1 <- ClassifyResult("Example", "Differential Expression", "t-test",
                            LETTERS[1:20], LETTERS[10:1], 100,
                            list(1:100, c(1:9, 11:101)), list(sample(10, 10), sample(10, 10)),
                            list(function(oracle){}), predicted, actual,
                            list("permuteFold", 2, 2))
  
  predicted[[1]][c(2, 6), "Healthy"] <- c(0.40, 0.60)
  predicted[[1]][c(2, 6), "Cancer"] <- c(0.60, 0.40)
  result2 <- ClassifyResult("Example", "Differential Variability", "Bartlett Test",
                            LETTERS[1:20], LETTERS[10:1], 100, list(1:100, c(1:5, 11:105)),
                            list(sample(10, 10), sample(10, 10)), list(function(oracle){}),
                            predicted, actual, validation = list("permuteFold", 2, 2))
  ROCplot(list(result1, result2), comparisonVariable = "classificationName",
          plotTitle = "Cancer ROC")
}