\name{getLocationsAndScales}
\alias{getLocationsAndScales}
\alias{getLocationsAndScales,matrix-method}
\alias{getLocationsAndScales,ExpressionSet-method}
\title{Calculate Location and Scale}
\description{Calculates the location and scale for each feature.
}
\usage{
  \S4method{getLocationsAndScales}{matrix}(expression, ...)
  \S4method{getLocationsAndScales}{ExpressionSet}(expression, location = c("mean", "median"),
                                   scale = c("SD", "MAD", "Qn"))
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    data. For a matrix, the rows are features, and the columns
                    are samples.}
  \item{...}{Unused variables from the \code{\link{matrix}} method passed to the
             \code{\link{ExpressionSet}} method.}
  \item{location}{The location to be calculated.}
  \item{scale}{The scale to be calculated.}
}
\details{
  Location can be either "mean" or "median". Scale can be standard deviation, median absolute deviation, or \eqn{Q_n}{Qn}.
}
\value{
  A \code{\link{list}} of length 2. The first element contains the location for
  every feature. The second element contains the scale for every feature.
}
\author{Dario Strbenac}
\references{
  Qn: \url{http://www.tandfonline.com/doi/pdf/10.1080/01621459.1993.10476408}
}
\examples{
  genesMatrix <- matrix(rnorm(1000, 8, 4), ncol = 10)
  getLocationsAndScales(genesMatrix, "median", "MAD")
}
