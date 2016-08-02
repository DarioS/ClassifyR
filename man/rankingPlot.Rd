\name{rankingPlot}
\alias{rankingPlot}
\alias{rankingPlot,list-method}
\title{Plot Pair-wise Overlap of Ranked Features}
\description{Pair-wise overlaps can be done for two types of analyses. Firstly, each cross-validation iteration
can be considered within a single classification. This explores the feature ranking stability. Secondly, the
overlap may be considered between different classification results. This approach compares the feature ranking
commonality between different methods. Two types of commonality are possible to analyse. One summary is
the average pair-wise overlap between a level of the comparison factor and the other summary is the pair-wise
overlap of each level of the comparison factor that is not the reference level against the reference level.
The overlaps are converted to percentages and plotted as lineplots.
       }
\usage{
  \S4method{rankingPlot}{list}(results, topRanked = seq(10, 100, 10),
                            comparison = c("within", "classificationName", "validation", "datasetName", "selectionName"),
                            referenceLevel = NULL,
                            lineColourVariable = c("validation", "datasetName", "classificationName",
                                                   "selectionName", "None"),
                            lineColours = NULL, lineWidth = 1,
                            pointTypeVariable = c("datasetName", "classificationName", "validation",
                                                  "selectionName", "None"),
                            pointSize = 2, legendLinesPointsSize = 1,
                            rowVariable = c("None", "datasetName", "classificationName", "validation", "selectionName"),
                            columnVariable = c("classificationName", "datasetName", "validation", "selectionName", "None"),
                            yMax = 100, fontSizes = c(24, 16, 12, 12, 12, 16),
                            title = if(comparison[1] == "within") "Feature Ranking Stability" else "Feature Ranking Commonality",
                            xLabelPositions = seq(10, 100, 10),
                            yLabel = if(is.null(referenceLevel)) "Average Common Features (\%)" else paste("Average Common Features with", referenceLevel, "(\%)"),
                            margin = grid::unit(c(0, 0, 0, 0), "lines"),
                            showLegend = TRUE, plot = TRUE, parallelParams = bpparam())
}
\arguments{
  \item{results}{A list of \code{\link{ClassifyResult}} or \code{\link{SelectResult}} objects.}
  \item{topRanked}{A sequence of thresholds of number of the best features to use for overlapping.}
  \item{comparison}{The aspect of the experimental design to compare. See \code{Details} section for a
                    detailed description.}
  \item{referenceLevel}{The level of the comparison factor to use as the reference to compare each
                        non-reference level to. If \code{NULL}, then each level has the
                        average pairwise overlap calculated to all other levels.}                    
  \item{lineColourVariable}{The slot name that different levels of are plotted as
                            different line colours.}
  \item{lineColours}{A vector of colours for different levels of the line colouring parameter. If \code{NULL},
                     a default palette is used.}       
  \item{lineWidth}{A single number controlling the thickness of lines drawn.}
  \item{pointTypeVariable}{The slot name that different levels of are plotted as
                            different point shapes on the lines.}
  \item{pointSize}{A single number specifying the diameter of points drawn.}
  \item{legendLinesPointsSize}{A single number specifying the size of the lines and points in the legend,
                               if a legend is drawn.}
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
  \item{margin}{The margin to have around the plot.}
  \item{showLegend}{If \code{TRUE}, a legend is plotted next to the plot. If FALSE, it is hidden.}  
  \item{plot}{Logical. If \code{TRUE}, a plot is produced on the current graphics device.}
  \item{parallelParams}{An object of class \code{\link{MulticoreParam}} or \code{\link{SnowParam}}.}
}
\details{
  Possible values for characteristics are \code{"datasetName"}, \code{"classificationName"},
  \code{"selectionName"}, and \code{"validation"}. If \code{"None"}, then that graphical element is not used.
  
  If \code{comparison} is \code{"within"}, then the feature rankings are compared within a particular
  analysis. The result will inform how stable the feature rankings are between different iterations of cross-validation for a particular analysis. If \code{comparison} is \code{"classificationName"}, then the feature
  rankings are compared across different classification algorithm types, for each level of \code{"datasetName"},
  \code{"selectionName"} and \code{"validation"}. The result will inform how stable the feature rankings
  are between different classification algorithms, for every cross-validation scheme, selction algorithm and
  dataset. If \code{comparison} is \code{"selectionName"}, then the feature rankings are compared across different
  feature selection algorithms, for each level of \code{"datasetName"}, \code{"classificationName"} and
  \code{"validation"}. The result will inform how stable the feature rankings are between feature selection
  classification algorithms, for every dataset, classification algorithm, and cross-validation scheme.
  If \code{comparison} is \code{"validation"}, then the feature rankings are compared across different
  cross-validation schemes, for each level of \code{"classificationName"}, \code{"selectionName"} and
  \code{"datasetName"}. The result will inform how stable the feature rankings are between different
  cross-validation schemes, for every selection algorithm, classification algorithm and every dataset.
  If \code{comparison} is \code{"datasetName"}, then the feature rankings are compared across different datasets,
  for each level of \code{"classificationName"}, \code{"selectionName"} and \code{"validation"}.
  The result will inform how stable the feature rankings are between different datasets, for every
  classification algorithm and every dataset. This could be used to consider if different experimental
  studies have a highly overlapping feature ranking pattern.
  
  Calculating all pair-wise set overlaps for a large cross-validation result can be time-consuming.
  This stage can be done on multiple CPUs by providing the relevant options to \code{parallelParams}.
}
\value{
  An object of class \code{ggplot} and a plot on the current graphics device, if \code{plot} is \code{TRUE}.
}
\author{Dario Strbenac}

\examples{
  predicted <- data.frame(sample = sample(10, 100, replace = TRUE),
                          label = rep(c("Healthy", "Cancer"), each = 50))
  actual <- factor(rep(c("Healthy", "Cancer"), each = 5))
  rankList <- list(list(1:100, c(5:1, 6:100)), list(c(1:9, 11:101), c(1:50, 60:51, 61:100)))
  result1 <- ClassifyResult("Example", "Differential Expression", "Example Selection", LETTERS[1:10], LETTERS[10:1],
                            rankList,
                            list(list(rankList[[1]][[1]][1:15], rankList[[1]][[2]][1:15]),
                                 list(rankList[[2]][[1]][1:10], rankList[[2]][[2]][1:10])),
                            list(predicted), actual, list("resampleFold", 2, 2))
  
  predicted[, "label"] <- sample(predicted[, "label"])
  rankList <- list(list(1:100, c(sample(20), 21:100)), list(c(1:9, 11:101), c(1:50, 60:51, 61:100)))
  result2 <- ClassifyResult("Example", "Differential Variability", "Example Selection", LETTERS[1:10], LETTERS[10:1],
                            rankList,
                            list(list(rankList[[1]][[1]][1:15], rankList[[1]][[2]][1:15]),
                                 list(rankList[[2]][[1]][1:10], rankList[[2]][[2]][1:10])),
                            list(predicted), actual, validation = list("resampleFold", 2, 2))
                            
  rankingPlot(list(result1, result2), pointTypeVariable = "classificationName")
  
  oneRanking <- c(10, 8, 1, 2, 3, 4, 7, 9, 5, 6)
  otherRanking <- c(8, 2, 3, 4, 1, 10, 6, 9, 7, 5)
  oneResult <- SelectResult("Example", "One Method", list(oneRanking), list(oneRanking[1:5]))
  otherResult <- SelectResult("Example", "Another Method", list(otherRanking), list(otherRanking[1:2]))
  
  rankingPlot(list(oneResult, otherResult), comparison = "selectionName",
           referenceLevel = "One Method", topRanked = seq(2, 8, 2),
           lineColourVariable = "selectionName", columnVariable = "None",
           pointTypeVariable = "None", xLabelPositions = 1:10)
}