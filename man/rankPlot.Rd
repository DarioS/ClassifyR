\name{rankPlot}
\alias{rankPlot}
\alias{rankPlot,list-method}
\title{Plot Pair-wise Overlap of Ranked Features}
\description{The average pair-wise overlap is computed for every pair of cross-validations. The overlap is converted to a percentage and plotted as lineplots.
       }
\usage{
  \S4method{rankPlot}{list}(results, topRanked = seq(10, 100, 10),
                            comparison = c("within", "classificationName", "validation", "datasetName"),
                            lineColourVariable = c("validation", "datasetName", "classificationName", "None"),
                            lineColours = NULL, lineWidth = 1,
                            pointTypeVariable = c("datasetName", "classificationName", "validation", "None"),
                            pointSize = 2, legendLinesPointsSize = 1,
                            rowVariable = c("None", "datasetName", "classificationName", "validation"),
                            columnVariable = c("classificationName", "datasetName", "validation", "None"),
                            yMax = 100, fontSizes = c(24, 16, 12, 12, 12, 16), title = "Feature Ranking Stability",
                            xLabelPositions = seq(10, 100, 10), yLabel = "Average Pairwise Common Features (\%)",
                            plot = TRUE, parallelParams = bpparam())
}
\arguments{
  \item{results}{A list of \code{\link{ClassifyResult}} objects.}
  \item{topRanked}{A sequence of thresholds of number of the best features to use for overlapping.}
  \item{comparison}{The aspect of the experimental design to compare. See \code{Details} section for a
                    detailed description.}  
  \item{lineColourVariable}{The slot name that different levels of are plotted as
                            different line colours.}
  \item{lineColours}{A vector of colours for different levels of the line colouring parameter. If \code{NULL},
                     a default palette is used.}       
  \item{lineWidth}{A single number controlling the thickness of lines drawn.}
  \item{pointTypeVariable}{The slot name that different levels of are plotted as
                            different point shapes on the lines.}
  \item{pointSize}{A single number diameter of points drawn.}
  \item{legendLinesPointsSize}{A single number specifying the size of the lines and points in the legend, if a legend is drawn.}
  \item{rowVariable}{The slot name that different levels of are plotted as separate rows of lineplots.}     
  \item{columnVariable}{The slot name that different levels of are plotted as separate columns of lineplots.}
  \item{yMax}{The maximum value of the percentage to plot.}
  \item{fontSizes}{A vector of length 6. The first number is the size of the title.
                   The second number is the size of the axes titles. The third number is
                   the size of the axes values. The fourth number is the size of the
                   legends' titles. The fifth number is the font size of the legend labels.
                   The sixth number is the font size of the titles of grouped plots, if any
                   are produced. In other words, when \code{rowVariable} or
                   \code{columnVariable} are not \code{NULL}.}
  \item{title}{An overall title for the plot.}
  \item{xLabelPositions}{Locations where to put labels on the x-axis.}
  \item{yLabel}{Label to be used for the y-axis of overlap percentages.}
  \item{plot}{Logical. IF TRUE, a plot is produced on the current graphics device.}
  \item{parallelParams}{An object of class \code{\link{MulticoreParam}} or \code{\link{SnowParam}}.}
}
\details{
  Possible values for slot names are \code{"datasetName"}, \code{"classificationName"}, and
  \code{"validation"}. If \code{"None"}, then that graphic element is not used.
  
  If \code{comparison} is \code{"within"}, then the feature rankings are compared within a particular
  analysis. The result will inform how stable the feature rankings are between different iterations of a
  particular analysis. If \code{comparison} is \code{"classificationName"}, then the feature
  rankings are compared across different classification algorithm types, for each level of \code{"datasetName"} and
  \code{"validation"}. The result will inform how stable the feature rankings are between different
  classification algorithms, for every cross-validation scheme and every dataset. If \code{comparison} is
  \code{"validation"}, then the feature rankings are compared across different cross-validation schemes, for
  each level of \code{"classificationName"} and \code{"datasetName"}. The result will inform how stable the
  feature rankings are between different cross-validation schemes, for every classification algorithm and
  every dataset. If \code{comparison} is \code{"datasetName"}, then the feature rankings are compared across
  different datasets, for each level of \code{"classificationName"} and \code{"validation"}.
  The result will inform how stable the feature rankings are between different datasets, for
  every classification algorithm and every dataset. This could be used to consider if different studies have
  a highly overlapping feature ranking pattern.
  
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
  result1 <- ClassifyResult("Example", "Differential Expression", LETTERS[1:10], LETTERS[10:1], list(1:100, c(1:9, 11:101)), list(sample(10, 10), sample(10, 10)),
                            list(predicted), actual, list("fold", 100, 5))
  predicted[, "label"] <- sample(predicted[, "label"])
  result2 <- ClassifyResult("Example", "Differential Variability", LETTERS[1:10], LETTERS[10:1], list(1:100, c(1:5, 11:105)), list(sample(10, 10), sample(10, 10)),
                            list(predicted), actual, validation = list("leave", 1))
  # rankPlot(list(result1, result2), pointTypeVariable = "classificationName") # Wait for namespace problems to be fixed.
}