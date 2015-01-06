\name{plotFeatureClasses}
\alias{plotFeatureClasses}
\alias{plotFeatureClasses,matrix-method}
\alias{plotFeatureClasses,ExpressionSet-method}
\title{Plot Density and Scatterplot for Genes By Class}
\description{Allows the visualisation of genes which were selected by a feature
selection method.
       }
\usage{
  \S4method{plotFeatureClasses}{matrix}(expression, classes, ...)
  \S4method{plotFeatureClasses}{ExpressionSet}(expression, rows, plot = c("both", "density", "stripchart"),
                               xAxisLabel = expression(log[2](expression)), expressionLimits = c(2, 16),
                               yAxisLabels = c("Density", "Classes"), showXtickLabels = TRUE,
                               showYtickLabels = TRUE, fontSizes = c(24, 16, 12, 12, 12),
                               colours = c("blue", "red"))
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    the training data. For a matrix, the rows are features, and the columns
                    are samples.}
  \item{classes}{A vector of class labels.}
  \item{...}{Unused variables from the \code{\link{matrix}} method passed to the
             \code{\link{ExpressionSet}} method.}
  \item{rows}{A vector specifying which rows of the matrix to plot.}             
  \item{plot}{Which plots to draw. Can draw either a density plot, stripchart, or both.}
  \item{xAxisLabel}{The axis label for the expression axis.}
  \item{yAxisLabels}{A character vector of length 2. The first value is the y-axis label
                     for the density plot. The second value is the y-axis labels for the
                     stripchart. Provide both labels, even if only plotting one kind of plot.}
  \item{expressionLimits}{The minimum and maximum expression values to plot. Set to NULL to use range of data.}
  \item{showXtickLabels}{Logical. IF FALSE, the x-axis labels are hidden.}
  \item{showYtickLabels}{Logical. IF FALSE, the y-axis labels are hidden.}  
  \item{fontSizes}{A vector of length 5. The first number is the size of the title.
                   The second number is the size of the axes titles. The third number is
                   the size of the axes values. The fourth number is the size of the
                   legends' titles. The fifth number is the font size of the legend labels.}  
  \item{colours}{The colours to plot data of each class in.}
}
\value{
  Plots.
}
\author{Dario Strbenac}

\examples{
  # First 25 samples are mixtures of two normals. Last 25 samples are one normal.
  genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(50, 5, 1), rnorm(50, 15, 1)))
  genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn) rnorm(100, 9, 3)))
  classes <- factor(rep(c("Poor", "Good"), each = 25), levels = c("Good", "Poor"))
  chosen <- 1:5 # First five genes in the data were chosen.
                                     
  # plotFeatureClasses(genesMatrix, classes, chosen, expressionLimits = NULL) Waiting for gridExtra fix.
}