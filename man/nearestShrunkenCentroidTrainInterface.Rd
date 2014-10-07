\name{nearestShrunkenCentroidTrainInterface}
\alias{nearestShrunkenCentroidTrainInterface}
\alias{nearestShrunkenCentroidTrainInterface,matrix-method}
\alias{nearestShrunkenCentroidTrainInterface,ExpressionSet-method}

\title{Interface for \code{pamr.train} Function from \code{pamr} CRAN Package}
\description{Restructures variables from ClassifyR framework to be compatible with \code{\link[pamr]{pamr.train}} definition.}
\usage{
  \S4method{nearestShrunkenCentroidTrainInterface}{matrix}(expression, classes, ...)
  \S4method{nearestShrunkenCentroidTrainInterface}{ExpressionSet}(expression, ..., verbose = 3)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    the training data. For a matrix, the rows are features, and the columns
                    are samples.}             
  \item{classes}{A vector of class labels.}                               
  \item{...}{Extra arguments passed to \code{\link[pamr]{pamr.train}}.}            
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  This function is an interface between the ClassifyR framework and \code{\link[pamr]{pamr.train}}.
}
\value{
  A list with elements as described in \code{\link[pamr]{pamr.train}}.
}
\author{Dario Strbenac}
\seealso{\code{\link[pamr]{pamr.train}} for the function that was interfaced to.}
\examples{
  if(require(pamr))
  {
    # Samples in one class with differential expression to other class.
    genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(100, 9, 1)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn)
                                 c(rnorm(75, 9, 1), rnorm(25, 14, 1))))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    
    nearestShrunkenCentroidTrainInterface(genesMatrix, classes)
  }
}
