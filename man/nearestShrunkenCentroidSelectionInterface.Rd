\name{nearestShrunkenCentroidSelectionInterface}
\alias{nearestShrunkenCentroidSelectionInterface}
\alias{nearestShrunkenCentroidSelectionInterface,matrix-method}
\alias{nearestShrunkenCentroidSelectionInterface,ExpressionSet-method}

\title{Interface for \code{pamr.listgenes} Function from \code{pamr} CRAN Package}
\description{Restructures variables from ClassifyR framework to be compatible with \code{\link[pamr]{pamr.listgenes}} definition.}
\usage{
  \S4method{nearestShrunkenCentroidSelectionInterface}{matrix}(expression, classes, ...)
  \S4method{nearestShrunkenCentroidSelectionInterface}{ExpressionSet}(expression, trained, ..., verbose = 3)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    the training data. For a matrix, the rows are features, and the columns
                    are samples.}             
  \item{classes}{A vector of class labels.} 
  \item{trained}{The output of \code{\link{nearestShrunkenCentroidTrainInterface}},
                      which is identical to the output of \code{\link[pamr]{pamr.listgenes}}.}                  
  \item{...}{Extra arguments passed to \code{\link[pamr]{pamr.listgenes}}}            
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  This function is an interface between the ClassifyR framework and \code{\link[pamr]{pamr.listgenes}}.
  
  The set of features chosen is the obtained by considering the range of thresholds provided to \code{\link{nearestShrunkenCentroidTrainInterface}} and using the threshold that obtains the lowest cross-validation error rate on the training set.
}
\value{
  A vector of numeric row indexes, of the features that were chosen after applying the threshold.
}
\author{Dario Strbenac}
\seealso{\code{\link[pamr]{pamr.listgenes}} for the function that was interfaced to.}
\examples{
  if(require(pamr))
  {
    # Samples in one class with differential expression to other class.
    genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(100, 9, 1)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn)
                                 c(rnorm(75, 9, 1), rnorm(25, 14, 1))))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    
    trained <- nearestShrunkenCentroidTrainInterface(genesMatrix, classes)
    nearestShrunkenCentroidSelectionInterface(genesMatrix, classes, trained)
  }
}
