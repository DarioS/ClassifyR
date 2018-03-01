\name{getLocationsAndScales}
\alias{getLocationsAndScales}
\alias{getLocationsAndScales,matrix-method}
\alias{getLocationsAndScales,DataFrame-method}
\alias{getLocationsAndScales,MultiAssayExperiment-method}
\title{Calculate Location and Scale}
\description{Calculates the location and scale for each feature.
}
\usage{
  \S4method{getLocationsAndScales}{matrix}(measurements, ...)
  \S4method{getLocationsAndScales}{DataFrame}(measurements, location = c("mean", "median"),
                      scale = c("SD", "MAD", "Qn"))
  \S4method{getLocationsAndScales}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
                      \code{\link{MultiAssayExperiment}} containing the data.
                      For a \code{matrix}, the rows are features, and the columns are samples.}
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 numeric variables from the clinical data table will be used.}
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method.}  
  \item{location}{The type of location to be calculated.}
  \item{scale}{The type of scale to be calculated.}
}
\details{
  \code{"SD"} is used to represent standard deviation and \code{"MAD"} is used to represent
  median absolute deviation.
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
  distributionInfo <- getLocationsAndScales(genesMatrix, "median", "MAD")
  
  mean(distributionInfo[["median"]]) # Typical median.
  mean(distributionInfo[["MAD"]]) # Typical MAD.
}
