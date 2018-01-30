\name{randomForestInterface}
\alias{randomForestInterface}
\alias{randomForestInterface,matrix-method}
\alias{randomForestInterface,DataFrame-method}
\alias{randomForestInterface,MultiAssayExperiment-method}
\title{
  An Interface for randomForest Package's randomForest Function
}
\description{
  A random forest classifier builds multiple decision trees and uses the predictions of the trees
  to determine a single prediction for each test sample.
}
\usage{
  \S4method{randomForestInterface}{matrix}(measurements, classes, test, ...)
  \S4method{randomForestInterface}{DataFrame}(measurements, classes, test, ..., verbose = 3)
  \S4method{randomForestInterface}{MultiAssayExperiment}(measurements, targets = names(measurements),
                                                         test, ...)
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
                     \code{\link{MultiAssayExperiment}} containing the training data.
                     For a \code{matrix}, the rows are features, and the columns are samples.}
  \item{classes}{Either a vector of class labels of class \code{\link{factor}} of the same length
                 as the number of samples in \code{measurements} or if the measurements are
                 of class \code{DataFrame} a character vector of length 1 containing the
                 column name in \code{measurement} is also permitted. Not used if \code{measurements}
                 is a \code{MultiAssayExperiment} object.}
  \item{test}{An object of the same class as \code{measurements} with no samples in common with
              \code{measurements} and the same number of features as it.}                 
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 integer variables from the clinical data table will be used.}
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method (e.g. \code{verbose}) or
             options which are accepted by the \code{link[randomForest]{randomForest}} function.}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}             
}
\details{
  If \code{measurements} is an object of class \code{MultiAssayExperiment}, the factor of sample
  classes must be stored in the DataFrame accessible by the \code{colData} function with
  column name \code{"class"}.
}
\value{
  An object of type \code{randomForest}. The predictions of the test set samples are stored in
  the list element named \code{"predicted"} of the \code{"test"} element.
}
\author{Dario Strbenac}
\examples{
  if(require(randomForest))
  {
    # Genes 76 to 100 have differential expression.
    genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
                                      c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix))
    rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
    selected <- rownames(genesMatrix)[91:100]
    trainingSamples <- c(1:20, 26:45)
    testingSamples <- c(21:25, 46:50)
    
    classified <- randomForestInterface(genesMatrix[, trainingSamples],
                                        classes[trainingSamples],
                                        genesMatrix[, testingSamples], ntree = 10)
    classified[["test"]][["predicted"]]
  }
}
