\name{calcPerformance}
\alias{calcPerformance}
\alias{calcPerformance,ClassifyResult-method}
\title{Add Performance Calculations to a ClassifyResult object}
\description{Annotates the results of calling \code{\link{runTests}} with different kinds of performance measures.
}
\usage{
  \S4method{calcPerformance}{ClassifyResult}(result, performanceType, ...)
}
\arguments{
  \item{result}{An object of class \code{\link{ClassifyResult}}.}
  \item{performanceType}{Either \code{"balanced"} or one of the options provided by \code{\link[ROCR]{performance}}.}
  \item{...}{Further arguments that may be used by \code{\link[ROCR]{performance}}.}
}
\details{
  If \code{\link{runTests}} was run in resampling mode, one performance measure is produced
  for every resampling. If the leave-out mode was used, then the predictions are
  concatenated, and one performance measure is calcuated for all predictions.
  
  Because ROCR only provides calculations for two-class classification, this function
  is only suitable for two-class classification performance measures.
}
\value{
  An updated \code{\linkS4class{ClassifyResult}} object, with new information in the \code{performance} slot.
}
\examples{
  predictTable <- data.frame(sample = 1:5,
                             predicted = factor(sample(LETTERS[1:2], 50, replace = TRUE)))
  actual <- factor(sample(LETTERS[1:2], 50, replace = TRUE))                             
  result <- ClassifyResult("Example", "Differential Expression",
                           paste("A", 1:10, sep = ''), paste("Gene", 1:50, sep = ''),
                           list(1:100, 1:100), list(1:5, 6:15),
                           list(predictTable), actual, list("leave", 2))
  result <- calcPerformance(result, "balanced") 
  performance(result)
}
\author{Dario Strbenac}
