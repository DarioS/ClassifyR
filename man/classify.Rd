\name{classify}
\alias{classify}
\title{
Interface for PoiClaClu Package's Classify Function
}
\description{
Passes along all parameters except \code{verbose}, from the framework to \code{\link[PoiClaClu]{Classify}}.
}
\usage{
classify(...)
}
\arguments{
  \item{...}{All parameters that \code{\link[PoiClaClu]{Classify}} can accept and also \code{verbose}.}
}
\value{
  A result list, the same as is returned by \code{\link[PoiClaClu]{Classify}}.
}
\author{Dario Strbenac}
\examples{
  if(require(PoiClaClu))
  {
    readCounts <- CountDataSet(n = 100, p = 1000, 2, 5, 1)
    classify(readCounts[["x"]], readCounts[["y"]], readCounts[["xte"]], verbose = TRUE)
  }
}
