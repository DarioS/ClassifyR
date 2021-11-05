\name{subtractFromLocation}
\alias{subtractFromLocation}
\alias{subtractFromLocation,matrix-method}
\alias{subtractFromLocation,DataFrame-method}
\alias{subtractFromLocation,MultiAssayExperiment-method}
\title{Subtract Numeric Feature Measurements from a Location}
\description{For each numeric feature, calculates the location, and subtracts
             all measurements from that location.
             }
\usage{
  \S4method{subtractFromLocation}{matrix}(measurements, training, location = c("mean", "median"),
                     absolute = TRUE, transformName = "Location Subtraction", verbose = 3)
  \S4method{subtractFromLocation}{DataFrame}(measurements, training, location = c("mean", "median"),
                     absolute = TRUE, transformName = "Location Subtraction", verbose = 3)
  \S4method{subtractFromLocation}{MultiAssayExperiment}(measurements, training, targets = names(measurements),
                           location = c("mean", "median"), absolute = TRUE, transformName = "Location Subtraction", verbose = 3)
}
\arguments{
  \item{measurements}{A \code{\link{matrix}}, \code{\link{DataFrame}} or a
                      \code{\link{MultiAssayExperiment}} object containing the data.
                      For a matrix, the rows are for features and the columns are
                      for samples.}
  \item{training}{A vector specifying which samples are in the training set.}                
  \item{location}{Character. Either "mean" or "median".}
  \item{absolute}{Logical. Default: \code{TRUE}. If \code{TRUE}, then absolute values of the differences
                  are returned. Otherwise, they are signed.}
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 numeric variables from the clinical data table will be used.}                   \item{transformName}{Default: Location Subtraction. Useful for automated plot annotation by
                       plotting functions within this package..}
  \item{verbose}{Default: 3. A progress message is shown if this value is 3.}
}
\details{
  Only the samples specified by \code{training} are used in the calculation of the location.
  To use all samples for calculation of the location, simply provide indices of all the samples.
}
\value{
  The same class of variable as the input variable \code{measurements} is, with the numeric features
  subtracted from the calculated location.
}
\author{Dario Strbenac}

\examples{
  aMatrix <- matrix(1:100, ncol = 10)
  subtractFromLocation(aMatrix, training = 1:5, "median")
}