\name{calcError}
\alias{calcError}
\alias{calcError,ClassifyResult-method}
\title{Add Error Calculations to a ClassifyResult object}
\description{Annotates the results of calling \code{\link{runTests}} with different kinds of error measures.
}
\usage{
  \S4method{calcError}{ClassifyResult}(result, errorType, ...)
}
\arguments{
  \item{result}{An object of class \code{\link{ClassifyResult}}.}
  \item{errorType}{Either "balanced" or one of the options provided by \code{\link[ROCR]{performance}}.}
  \item{...}{Further arguments that may be used by \code{\link[ROCR]{performance}}.}
}
\details{
  If \code{\link{runTests}} was run in resampling mode, one error measure is produced
  for every resampling. If the leave-out mode was used, then the predictions are
  concatenated, and one error measure is calcuated for all predictions.
  
  Because ROCR only provides calculations for two-class classification, this function
  is only suitable for two-class classification error measures.
}
\value{
  An updated \code{\linkS4class{ClassifyResult}} object, with new information in the \code{errors} slot.
}
\examples{
  predictTable <- data.frame(sample = 1:5,
                             predicted = factor(sample(LETTERS[1:2], 50, replace = TRUE)))
  actual <- factor(sample(LETTERS[1:2], 50, replace = TRUE))                             
  result <- ClassifyResult(paste("A", 1:10, sep = ''), paste("Gene", 1:50, sep = ''),
                           list(1:5, 6:15), list(predictTable), actual, list("leave", 2))
  result <- calcError(result, "balanced") 
  errors(result)
}
\author{Dario Strbenac}
