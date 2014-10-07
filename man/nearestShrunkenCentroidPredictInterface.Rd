\name{nearestShrunkenCentroidPredictInterface}
\alias{nearestShrunkenCentroidPredictInterface}
\alias{nearestShrunkenCentroidPredictInterface,pamrtrained,matrix-method}
\alias{nearestShrunkenCentroidPredictInterface,pamrtrained,ExpressionSet-method}

\title{Interface for \code{pamr.predict} Function from \code{pamr} CRAN Package}
\description{Restructures variables from ClassifyR framework to be compatible with \code{\link[pamr]{pamr.predict}} definition.}
\usage{
  \S4method{nearestShrunkenCentroidPredictInterface}{pamrtrained,matrix}(trained, test, ...)
  \S4method{nearestShrunkenCentroidPredictInterface}{pamrtrained,ExpressionSet}(trained, test, ..., verbose = 3)
}
\arguments{
  \item{trained}{An object of class \code{pamrtrained}.}             
  \item{test}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
              the test data. For a matrix, the rows are features, and the columns are samples.}            
  \item{...}{For the function with \code{\link{matrix}} dispatch, arguments passed to the function
             with \code{\link{ExpressionSet}} dispatch. For the function with \code{\link{ExpressionSet}}
             dispatch, arguemnts passed to \code{\link[pamr]{pamr.predict}}.} 
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  This function is an interface between the ClassifyR framework and \code{\link[pamr]{pamr.predict}}.
}
\value{
  A factor of predicted classes for the test data.
}
\author{Dario Strbenac}
\seealso{\code{\link[pamr]{pamr.predict}} for the function that was interfaced to.}
\examples{
  if(require(pamr))
  {
    # Samples in one class with differential expression to other class.
    genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(100, 9, 1)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn)
                                 c(rnorm(75, 9, 1), rnorm(25, 14, 1))))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    
    fit <- nearestShrunkenCentroidTrainInterface(genesMatrix[, c(1:20, 26:45)], classes[c(1:20, 26:45)])
    nearestShrunkenCentroidPredictInterface(fit, genesMatrix[, c(21:25, 46:50)])
  }
}
