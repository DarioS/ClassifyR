\name{selectionPlot}
\alias{selectionPlot}
\alias{selectionPlot,list-method}
\title{Plot Pair-wise Overlap or Selection Size Distribution of Selected Features}
\description{Pair-wise overlaps can be done for two types of analyses. Firstly, each cross-validation iteration
can be considered within a single classification. This explores the feature selection stability. Secondly, the
overlap may be considered between different classification results. This approach compares the feature selection
commonality between different selection methods. Two types of commonality are possible to analyse. One summary is
the average pair-wise overlap between all levels of the comparison factor and the other summary is the pair-wise
overlap of each level of the comparison factor that is not the reference level against the reference level.
The overlaps are converted to percentages and plotted as lineplots.

Additionally, a heatmap of selection size frequencies can be made by specifying size as the
comparison to make.
       }
\usage{
  \S4method{selectionPlot}{list}(results,
              comparison = "within", referenceLevel = NULL,
              characteristicsList = list(x = "Classifier Name"), coloursList = list(),
              orderingList = list(), binsList = list(), yMax = 100, fontSizes = c(24, 16, 12, 16),
              title = if(comparison[1] == "within") "Feature Selection Stability"
                      else if(comparison == "size") "Feature Selection Size" else
                      "Feature Selection Commonality",
              yLabel = if(is.null(referenceLevel) && comparison != "size") "Common Features (\%)"
                       else if(comparison == "size") "Set Size" else
                       paste("Common Features with", referenceLevel, "(\%)"),
              margin = grid::unit(c(1, 1, 1, 1), "lines"), rotate90 = FALSE,
              showLegend = TRUE, plot = TRUE, parallelParams = bpparam())
}
\arguments{
  \item{results}{A list of \code{\link{ClassifyResult}} objects.}
  \item{comparison}{Default: within. The aspect of the experimental design to compare. Can be any
                    characteristic that all results share or either one of the special values \code{"within"}
                    to compare between all pairwise iterations of cross-validation. or \code{"size"},
                    to draw a bar chart of the frequency of selected set sizes.}
  \item{referenceLevel}{The level of the comparison factor to use as the reference to compare each
                        non-reference level to. If \code{NULL}, then each level has the
                        average pairwise overlap calculated to all other levels.}                    
  \item{characteristicsList}{A named list of characteristics. Each element's name must be one of
  \code{"x"}, \code{"row"}, \code{"column"}, \code{fillColour}, or \code{fillLine}. The value of each
  element must be a characteristic name, as stored in the \code{"characteristic"} column of
  the results' characteristics table. Only \code{"x"} is mandatory.}
  \item{coloursList}{A named list of plot aspects and colours for the aspects. No elements are
  mandatory. If specified, each list element's name must be either \code{"fillColours"} or
  \code{"lineColours"}. If a characteristic is associated to fill or line by
  \code{characteristicsList} but this list is empty, a palette of colours will be
  automaticaly chosen.}
    \item{orderingList}{An optional named list. Any of the variables specified to \code{characteristicsList} can be the name of an element of this list and the value of
    the element is the order in which the factors should be presented in, in case
    alphabetical sorting is undesirable.}
    \item{binsList}{Used only if \code{comparison} is \code{"size"}. A list with elements named \code{"setSizes"}
    and \code{"frequencies"} Both elements are mandatory. \code{"setSizes"} specifies the bin boundaries for
    bins of interest of feature selection sizes (e.g. 0, 10, 20, 30). \code{"frequencies"} specifies the bin
    boundaries for the relative frequency percentages to plot (e.g. 0, 20, 40, 60, 80, 100).}
  \item{yMax}{Used only if \code{comparison} is not \code{"size"}. The maximum value of the percentage
    overlap to plot.}
  \item{fontSizes}{A vector of length 4. The first number is the size of the title.
                   The second number is the size of the axes titles. The third number is
                   the size of the axes values. The fourth number is the font size of the
                   titles of grouped plots, if any are produced. In other words, when
                   \code{rowVariable} or \code{columnVariable} are not \code{NULL}.}
  \item{title}{An overall title for the plot. By default, specifies whether stability or commonality is shown.}
  \item{yLabel}{Label to be used for the y-axis of overlap percentages. By default, specifies whether
  stability or commonality is shown.}
  \item{margin}{The margin to have around the plot.}
  \item{rotate90}{Logical. If \code{TRUE}, the boxplot is horizontal.}
  \item{showLegend}{If \code{TRUE}, a legend is plotted next to the plot. If FALSE, it is hidden.}  
  \item{plot}{Logical. If \code{TRUE}, a plot is produced on the current graphics device.}
  \item{parallelParams}{An object of class \code{\link{MulticoreParam}} or \code{\link{SnowParam}}.}
}
\details{
  If \code{comparison} is \code{"within"}, then the feature selection overlaps are compared within a particular
  analysis. The result will inform how stable the selections are between different iterations of cross-validation
  for a particular analysis. Otherwise, the comparison is between different cross-validation runs, and this gives
  an indication about how common are the features being selected by different classifications.
  
  Calculating all pair-wise set overlaps can be time-consuming. This stage can be done on multiple CPUs by
  providing the relevant options to \code{parallelParams}. The percentage is calculated as the intersection
  of two sets of features divided by the union of the sets, multiplied by 100.
  
  For the feature selection size mode, \code{binsList} is used to create bins which include the lowest value for
  the first bin, and the highest value for the last bin using \code{\link{cut}}.
}
\value{
  An object of class \code{ggplot} and a plot on the current graphics device, if \code{plot} is \code{TRUE}.
}
\author{Dario Strbenac}

\examples{
  predicted <- data.frame(sample = sample(10, 100, replace = TRUE),
                          class = rep(c("Healthy", "Cancer"), each = 50))
  actual <- factor(rep(c("Healthy", "Cancer"), each = 5))
  allFeatures <- sapply(1:100, function(index) paste(sample(LETTERS, 3), collapse = ''))
  rankList <- list(allFeatures[1:100], allFeatures[c(5:1, 6:100)],
                   allFeatures[c(1:9, 11, 10, 12:100)], allFeatures[c(1:50, 60:51, 61:100)])
  result1 <- ClassifyResult(DataFrame(characteristic = c("Data Set", "Selection Name", "Classifier Name",
                                                         "Cross-validations"),
                            value = c("Melanoma", "t-test", "Random Forest", "2 Permutations, 2 Folds")),
                            LETTERS[1:10], allFeatures, rankList,
                            list(rankList[[1]][1:15], rankList[[2]][1:15],
                                 rankList[[3]][1:10], rankList[[4]][1:10]),
                            list(function(oracle){}), NULL,
                            predicted, actual)
  
  predicted[, "class"] <- sample(predicted[, "class"])
  rankList <- list(allFeatures[1:100], allFeatures[c(sample(20), 21:100)],
                   allFeatures[c(1:9, 11, 10, 12:100)], allFeatures[c(1:50, 60:51, 61:100)])
  result2 <- ClassifyResult(DataFrame(characteristic = c("Data Set", "Selection Name", "Classifier Name",
                                                         "Cross-validation"),
                            value = c("Melanoma", "t-test", "Diagonal LDA", "2 Permutations, 2 Folds")),
                            LETTERS[1:10], allFeatures, rankList,
                            list(rankList[[1]][1:15], rankList[[2]][1:25],
                                 rankList[[3]][1:10], rankList[[4]][1:10]),
                            list(function(oracle){}), NULL,
                            predicted, actual)
  cList <- list(x = "Classifier Name", fillColour = "Classifier Name")
  selectionPlot(list(result1, result2), characteristicsList = cList)
  
  cList <- list(x = "Classifier Name", fillColour = "size")
  selectionPlot(list(result1, result2), comparison = "size",
                characteristicsList = cList,
                binsList = list(frequencies = seq(0, 100, 10), setSizes = seq(0, 25, 5))
                )
                
}