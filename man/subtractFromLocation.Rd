\name{subtractFromLocation}
\alias{subtractFromLocation}
\alias{subtractFromLocation,matrix-method}
\alias{subtractFromLocation,ExpressionSet-method}
\title{Subtract All Feature Measurements from Location}
\description{For each feature, calculates the location, and subtracts all measurements from
that location.
       }
\usage{
  \S4method{subtractFromLocation}{matrix}(expression, ...)
  \S4method{subtractFromLocation}{ExpressionSet}(expression, training, location = c("mean", "median"),
                                  absolute = TRUE, verbose = 3)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    the data. For a matrix, the rows are features, and the columns
                    are samples.}
  \item{...}{Unused variables from the \code{\link{matrix}} method passed to the
             \code{\link{ExpressionSet}} method.}                    
  \item{training}{A vector specifying which samples are in the training set.}                   
  \item{location}{Character. Either "mean" or "median".}
  \item{absolute}{If \code{TRUE}, then absolute values of the differences are returned.}
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 A higher number will produce more messages.}
}
\details{
  Only the samples specified by \code{training} are used in the calculation of the location.
  To use all samples for calculation of the location, simply provide indices of all the samples.
}
\value{
  An \code{\link{ExpressionSet}} of the same dimension that was input, with values
  subtracted from the location specified.
}
\author{Dario Strbenac}

\examples{
  subtractFromLocation(matrix(1:100, ncol = 10), training = 1:5, "median")
}