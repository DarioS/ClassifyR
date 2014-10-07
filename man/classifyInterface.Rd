\name{classifyInterface}
\alias{classifyInterface}
\title{
Interface for PoiClaClu Package's Classify Function
}
\description{
Passes along all parameters except \code{verbose}, from the framework to \code{\link[PoiClaClu]{Classify}}.
}
\usage{
classifyInterface(..., verbose = 3)
}
\arguments{
  \item{...}{All parameters that \code{\link[PoiClaClu]{Classify}} can accept and also \code{verbose}.}
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints a progress message if the value is 3.}
}
\value{
  A result list, the same as is returned by \code{\link[PoiClaClu]{Classify}}.
}
\author{Dario Strbenac}
\examples{
  if(require(PoiClaClu))
  {
    readCounts <- CountDataSet(n = 100, p = 1000, 2, 5, 1)
    classifyInterface(readCounts[["x"]], readCounts[["y"]], readCounts[["xte"]], verbose = TRUE)
  }
}
