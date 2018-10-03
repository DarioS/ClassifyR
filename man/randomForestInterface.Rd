\name{randomForestInterface}
\alias{randomForestInterface}
\alias{randomForestTrainInterface}
\alias{randomForestTrainInterface,matrix-method}
\alias{randomForestTrainInterface,DataFrame-method}
\alias{randomForestTrainInterface,MultiAssayExperiment-method}
\alias{randomForestPredictInterface}
\alias{randomForestPredictInterface,randomForest,matrix-method}
\alias{randomForestPredictInterface,randomForest,DataFrame-method}
\alias{randomForestPredictInterface,randomForest,MultiAssayExperiment-method}
\title{
  An Interface for randomForest Package's randomForest Function
}
\description{
  A random forest classifier builds multiple decision trees and uses the predictions of the trees
  to determine a single prediction for each test sample.
}
\usage{
  \S4method{randomForestTrainInterface}{matrix}(measurements, classes, ...)
  \S4method{randomForestTrainInterface}{DataFrame}(measurements, classes, ..., verbose = 3)
  \S4method{randomForestTrainInterface}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
  \S4method{randomForestPredictInterface}{randomForest,matrix}(forest, test, ...)
  \S4method{randomForestPredictInterface}{randomForest,DataFrame}(forest, test, ...,
                                          returnType = c("class", "score", "both"), verbose = 3)
  \S4method{randomForestPredictInterface}{randomForest,MultiAssayExperiment}(forest, test, targets = names(test), ...)
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
  \item{forest}{A trained random forest classifier, as created by
                \code{randomForestTrainInterface}, which has the same form as the
                output of \code{\link[randomForest]{randomForest}}.}
  \item{test}{An object of the same class as \code{measurements} with no samples in common with
              \code{measurements} and the same number of features as it.}                 
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 integer variables from the clinical data table will be used.}
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method (e.g. \code{verbose}) or
             options which are accepted by the \code{\link[randomForest]{randomForest}} or
             \code{\link[randomForest]{predict.randomForest}} functions.}
  \item{returnType}{Default: \code{"class"}. Either \code{"class"}, \code{"score"} or \code{"both"}.
                    Sets the return value from the prediction to either a vector of class labels,
                    score for a sample belonging to the second class, as determined by the
                    factor levels, or both labels and scores in a \code{data.frame}.}             
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}             
}
\details{
  If \code{measurements} is an object of class \code{MultiAssayExperiment}, the factor of sample
  classes must be stored in the DataFrame accessible by the \code{colData} function with
  column name \code{"class"}.
}
\value{
  For \code{randomForestTrainInterface}, the trained random forest. For \code{randomForestPredictInterface},
  either a factor vector of predicted classes, a vector of secores for
  the second class, or a table of both the class labels and second class scores,
  depending on the setting of \code{returnType}.
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
    colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix), sep = '')
    rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix), sep = '')
    trainingSamples <- c(1:20, 26:45)
    testingSamples <- c(21:25, 46:50)
    
    trained <- randomForestTrainInterface(genesMatrix[, trainingSamples],
                                          classes[trainingSamples])
    predicted <- randomForestPredictInterface(trained, genesMatrix[, testingSamples])
  }
}
