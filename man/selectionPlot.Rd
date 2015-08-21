\name{selectionPlot}
\alias{selectionPlot}
\alias{selectionPlot,list-method}
\title{Plot Pair-wise Overlap or Selection Size Distribution of Selected Features}
\description{Pair-wise overlaps can be done for two types of analyses. Firstly, each cross-validation iteration
can be considered within a single classification. This explores the feature selection stability. Secondly, the
overlap may be considered between different classification results. This approach compares the feature selection
commonality between different selection methods. Two types of commonality are possible to analyse. One summary is
the average pair-wise overlap between a level of the comparison factor and the other summary is the pair-wise
overlap of each level of the comparison factor that is not the reference level against the reference level.
The overlaps are converted to percentages and plotted as lineplots.

Additionally, a heatmap of selection size frequencies can be made.
       }
\usage{
  \S4method{selectionPlot}{list}(results,
                   comparison = c("within", "size", "classificationName", "validation", "datasetName", "selectionName"),
                   referenceLevel = NULL,
                   xVariable = c("classificationName", "datasetName", "validation", "selectionName"),
                   boxFillColouring = c("classificationName", "size", "datasetName", "validation",
                                        "selectionName", "None"),
                   boxFillColours = NULL,
                   boxFillBinBoundaries = NULL, setSizeBinBoundaries = NULL,
                   boxLineColouring = c("validation", "classificationName", "datasetName", "selectionName", "None"),
                   boxLineColours = NULL,
                   rowVariable = c("None", "validation", "datasetName", "classificationName", "selectionName"),
                   columnVariable = c("datasetName", "classificationName", "validation", "selectionName", "None"),
                   yMax = 100, fontSizes = c(24, 16, 12, 16),
                   title = if(comparison[1] == "within") "Feature Selection Stability" else if(comparison == "size") "Feature Selection Size" else "Feature Selection Commonality",
                   xLabel = "Analysis",
                   yLabel = if(is.null(referenceLevel) && comparison != "size") "Common Features (\%)" else if(comparison == "size") "Set Size" else paste("Common Features with", referenceLevel, "(\%)"),
                   margin = grid::unit(c(0, 1, 1, 0), "lines"), rotate90 = FALSE, plot = TRUE,
                   parallelParams = bpparam())
}
\arguments{
  \item{results}{A list of \code{\link{ClassifyResult}} or \code{\link{SelectResult}} objects.}
  \item{comparison}{The aspect of the experimental design to compare. See \code{Details} section for a
                    detailed description.}
  \item{referenceLevel}{The level of the comparison factor to use as the reference to compare each
                        non-reference level to. If \code{NULL}, then each level has the
                        average pairwise overlap calculated to all other levels.}                    
  \item{xVariable}{The factor to make separate boxes in the boxplot for.}
  \item{boxFillColouring}{A factor to colour the boxes by.}
  \item{boxFillColours}{A vector of colours, one for each level of \code{boxFillColouring}. If \code{NULL},
                     a default palette is used.}
  \item{boxFillBinBoundaries}{Used only if \code{comparison} is \code{"size"}. A vector of integers, specifying the bin
                     boundaries of percentages of size bins observed. e.g. 0, 10, 20, 30, 40, 50.}
  \item{setSizeBinBoundaries}{Used only if \code{comparison} is \code{"size"}. A vector of integers, specifying the bin
                     boundaries of set size bins. e.g. 50, 100, 150, 200, 250.}                         
  \item{boxLineColouring}{A factor to colour the box lines by.}
  \item{boxLineColours}{A vector of colours, one for each level of \code{boxLineColouring}. If \code{NULL},
                     a default palette is used.}
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
  \item{rotate90}{Logical. IF TRUE, the boxplot is horizontal.}
  \item{plot}{Logical. IF TRUE, a plot is produced on the current graphics device.}
  \item{parallelParams}{An object of class \code{\link{MulticoreParam}} or \code{\link{SnowParam}}.}
}
\details{
  Possible values for characteristics are \code{"datasetName"}, \code{"classificationName"}, \code{"size"},
  \code{"selectionName"}, and \code{"validation"}. If \code{"None"}, then that graphical element is not used.
  
  If \code{comparison} is \code{"within"}, then the feature selection overlaps are compared within a particular
  analysis. The result will inform how stable the selections are between different iterations of cross-validation
  for a particular analysis. If \code{comparison} is \code{"classificationName"}, then the feature
  selections are compared across different classification algorithm types, for each level of \code{"datasetName"},
  \code{"selectionName"} and \code{"validation"}. The result will inform how stable the feature selections
  are between different classification algorithms, for every cross-validation scheme, selction algorithm and
  dataset. If \code{comparison} is \code{"selectionName"}, then the feature selections are compared across
  different feature selection algorithms, for each level of \code{"datasetName"}, \code{"classificationName"} and
  \code{"validation"}. The result will inform how stable the feature selections are between feature selection
  algorithms, for every dataset, classification algorithm, and cross-validation scheme. If \code{comparison}
  is \code{"validation"}, then the feature selections are compared across different cross-validation schemes,
  for each level of \code{"classificationName"}, \code{"selectionName"} and \code{"datasetName"}. The result
  will inform how stable the feature selections are between different cross-validation schemes, for every
  selection algorithm, classification algorithm and every dataset. If \code{comparison} is
  \code{"datasetName"}, then the feature selections are compared across different datasets,
  for each level of \code{"classificationName"}, \code{"selectionName"}, and \code{"validation"}.
  The result will inform how stable the feature selections are between different datasets, for every
  classification algorithm and every dataset. This could be used to consider if different
  experimental studies have a highly overlapping feature selection pattern.
  
  Calculating all pair-wise set overlaps can be time-consuming. This stage can be done on multiple CPUs by
  providing the relevant options to \code{parallelParams}. The percentage is calculated as the intersection
  of two sets of features divided by the union of the sets, multiplied by 100.
  
  For the selection size mode, \code{boxFillBins} is used to create bins which include the lowest value for
  the first bin, and the highest value for the last bin using \code{\link{cut}}.
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
                            list(predicted), actual, list("fold", 2, 2))
  
  predicted[, "label"] <- sample(predicted[, "label"])
  rankList <- list(list(1:100, c(sample(20), 21:100)), list(c(1:9, 11:101), c(1:50, 60:51, 61:100)))
  result2 <- ClassifyResult("Example", "Differential Variability", "Example Selection", LETTERS[1:10], LETTERS[10:1],
                            rankList,
                            list(list(rankList[[1]][[1]][1:15], rankList[[1]][[2]][1:15]),
                                 list(rankList[[2]][[1]][1:10], rankList[[2]][[2]][1:10])),
                            list(predicted), actual, validation = list("fold", 2, 2))
                            
  selectionPlot(list(result1, result2), xVariable = "classificationName", xLabel = "Analysis", columnVariable = "None", rowVariable = "None", boxFillColouring = "classificationName")
  
  selectionPlot(list(result1, result2), comparison = "size", xVariable = "classificationName", xLabel = "Analysis", columnVariable = "None", rowVariable = "None", boxFillColouring = "size", boxFillBinBoundaries = seq(0, 100, 10),
  setSizeBinBoundaries = seq(0, 25, 5), boxLineColouring = "None")
  
  oneRanking <- c(10, 8, 1, 2, 3, 4, 7, 9, 5, 6)
  otherRanking <- c(8, 2, 3, 4, 1, 10, 6, 9, 7, 5)
  oneResult <- SelectResult("Example", "One Method", list(oneRanking), list(oneRanking[1:5]))
  otherResult <- SelectResult("Example", "Another Method", list(otherRanking), list(otherRanking[1:2]))
  
  selectionPlot(list(oneResult, otherResult), comparison = "selectionName", xVariable = "selectionName", xLabel = "Selection Method", columnVariable = "None", rowVariable = "None", boxFillColouring = "selectionName", boxLineColouring = "None")
}