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
  \S4method{plotFeatureClasses}{ExpressionSet}(expression, rows, assayType)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    the training data. For a matrix, the rows are features, and the columns
                    are samples.}
  \item{classes}{A vector of class labels.}
  \item{...}{Unused variables from the \code{\link{matrix}} method passed to the
             \code{\link{ExpressionSet}} method.}
  \item{rows}{A vector specifying which rows of the matrix to plot.}             
  \item{assayType}{Character. Either "array" or "sequencing".}
}
\value{
  Plots.
}
\author{Dario Strbenac}

\examples{
  # First 25 samples are mixtures of two normals. Last 25 samples are one normal.
  genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(50, 5, 1), rnorm(50, 15, 1)))
  genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn) rnorm(100, 9, 3)))
  classes <- factor(rep(c("Poor", "Good"), each = 25))
  chosen <- 1:5 # First five genes in the data were chosen.
                                     
  #plotFeatureClasses(genesMatrix, classes, chosen, "array")
}