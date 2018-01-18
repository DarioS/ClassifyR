\name{classifyInterface}
\alias{classifyInterface}
\alias{classifyInterface,matrix-method}
\alias{classifyInterface,DataFrame-method}
\alias{classifyInterface,MultiAssayExperiment-method}
\title{
  An Interface for PoiClaClu Package's Classify Function
}
\description{
  More details of Poisson LDA are available in the documentation of \code{\link[PoiClaClu]{Classify}}.
}
\usage{
  \S4method{classifyInterface}{matrix}(measurements, classes, test, ...)
  \S4method{classifyInterface}{DataFrame}(measurements, classes, test, ...)
  \S4method{classifyInterface}{MultiAssayExperiment}(measurements, test, targets = names(measurements), ...)
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
                     \code{\link{MultiAssayExperiment}} containing the training data.
                     For a \code{matrix}, the rows are features, and the columns are samples.
                     If of type \code{\link{DataFrame}}, the dataset is subset to only those
                     features of type \code{integer}.}
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
  \item{...}{All parameters that \code{\link[PoiClaClu]{Classify}} can accept.}
}
\value{
  A result list, the same as is returned by \code{\link[PoiClaClu]{Classify}}.
}
\author{Dario Strbenac}
\examples{
  if(require(PoiClaClu))
  {
    readCounts <- CountDataSet(n = 100, p = 1000, 2, 5, 0.1)
    # Rows are for features, columns are for samples.
    trainData <- t(readCounts[['x']])
    classes <- readCounts[['y']]
    testData <- t(readCounts[['xte']])
    classified <- classifyInterface(trainData, classes, testData, verbose = TRUE)
    
    setNames(table(readCounts[["yte"]] == classified[["ytehat"]]),
             c("Incorrect", "Correct"))
  }
}
