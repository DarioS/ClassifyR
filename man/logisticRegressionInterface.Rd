\name{logisticRegressionInterface}
\alias{logisticRegressionTrainInterface}
\alias{logisticRegressionTrainInterface,matrix-method}
\alias{logisticRegressionTrainInterface,DataFrame-method}
\alias{logisticRegressionTrainInterface,MultiAssayExperiment-method}
\alias{logisticRegressionPredictInterface}
\alias{logisticRegressionPredictInterface,mnlogit,matrix-method}
\alias{logisticRegressionPredictInterface,mnlogit,DataFrame-method}
\alias{logisticRegressionPredictInterface,mnlogit,MultiAssayExperiment-method}
\title{
  An Interface for mnlogit Package's mnlogit Function
}
\description{
  \code{logisticRegressionTrainInterface} generates a multinomial logistic regression model trained on some
  training data and \code{logisticRegressionPredictInterface} makes class predictions for samples in 
  the test dataset.
}
\usage{
  \S4method{logisticRegressionTrainInterface}{matrix}(measurements, classes, ...)
  \S4method{logisticRegressionTrainInterface}{DataFrame}(measurements, classes, ...)
  \S4method{logisticRegressionTrainInterface}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
  \S4method{logisticRegressionPredictInterface}{mnlogit,matrix}(model, test, ...)
  \S4method{logisticRegressionPredictInterface}{mnlogit,DataFrame}(model, test, ...)
  \S4method{logisticRegressionPredictInterface}{mnlogit,MultiAssayExperiment}(model, test, targets = names(test), ...)  
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
              \code{measurements} and the same number of features as it. Also, if a \code{DataFrame},
              the \code{class} column must be absent.}                 
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 integer variables from the clinical data table will be used.}
  \item{...}{Unused variables by the three top-level methods passed to the internal method
             which does the training or prediction (i.e. \code{verbose}).}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}  
}
\details{
  If \code{measurements} is an object of class \code{MultiAssayExperiment}, the factor of
  sample classes must be stored in the DataFrame accessible by the \code{colData} function
  with column name \code{"class"}.
  
  This wrapper works with individual-specific variables. If more a complex experimental design is utilised,
  such as a market research dataset with both individual-specific and alternative-specific variables, then
  this wrapper is not suitable to classify it.
}
\value{
  For \code{logisticRegressionTrainInterface}, a fitted multinomial logistic regression model.
  For \code{logisticRegressionPredictInterface}, a \code{factor} vector with class predictions for
  the samples in the test set.
}
\author{Dario Strbenac}
\examples{
  if(require(mnlogit))
  {
    variables <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
    trainSamples <- c(1:45, 51:95, 101:145)
    testSamples <- c(46:50, 96:100, 146:150)
    
    trained <- logisticRegressionTrainInterface(DataFrame(iris[, variables]),
                                                iris[, "Species"])
    predicted <- logisticRegressionPredictInterface(trained,
                                         DataFrame(iris[testSamples, variables]))
  }
}
