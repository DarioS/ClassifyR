\name{nearestShrunkenCentroidPredictInterface}
\alias{nearestShrunkenCentroidPredictInterface}
\alias{nearestShrunkenCentroidPredictInterface,pamrtrained,matrix-method}
\alias{nearestShrunkenCentroidPredictInterface,pamrtrained,DataFrame-method}
\alias{nearestShrunkenCentroidPredictInterface,pamrtrained,MultiAssayExperiment-method}

\title{Interface for \code{pamr.predict} Function from \code{pamr} CRAN Package}
\description{Restructures variables from ClassifyR framework to be compatible with \code{\link[pamr]{pamr.predict}} definition.}
\usage{
  \S4method{nearestShrunkenCentroidPredictInterface}{pamrtrained,matrix}(trained, test, ...)
  \S4method{nearestShrunkenCentroidPredictInterface}{pamrtrained,DataFrame}(trained, test, ...)  
  \S4method{nearestShrunkenCentroidPredictInterface}{pamrtrained,MultiAssayExperiment}(trained, test, targets = names(test), ...)
}
\arguments{
  \item{trained}{An object of class \code{pamrtrained}.}             
  \item{test}{An object of the same class as \code{measurements} with no samples in common with
              \code{measurements} used in the training stage and the same number of features as it.
              Also, if a \code{DataFrame}, the \code{class} column must be absent.}
  \item{targets}{If \code{test} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 numeric variables from the clinical data table will be used.}   
  \item{...}{Unused variables by the three top-level methods passed to the internal method
             which does the prediction or extra arguments passed to passed to \code{\link[pamr]{pamr.predict}}.} 
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
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
