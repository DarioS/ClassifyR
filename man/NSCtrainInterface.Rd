\name{NSCtrainInterface}
\alias{NSCtrainInterface}
\alias{NSCtrainInterface,matrix-method}
\alias{NSCtrainInterface,DataFrame-method}
\alias{NSCtrainInterface,MultiAssayExperiment-method}

\title{Interface for \code{pamr.train} Function from \code{pamr} CRAN Package}
\description{Restructures variables from ClassifyR framework to be compatible with
             \code{\link[pamr]{pamr.train}} definition.}
\usage{
  \S4method{NSCtrainInterface}{matrix}(measurements, classes, ...)
  \S4method{NSCtrainInterface}{DataFrame}(measurements, classes, ..., verbose = 3)
  \S4method{NSCtrainInterface}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
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
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 numeric variables from the clinical data table will be used.}                 
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method or extra arguments passed
             to \code{\link[pamr]{pamr.train}}.}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
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
    
    NSCtrainInterface(genesMatrix, classes)
  }
}
