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
  \item{performanceType}{Either \code{"balanced"}, \code{"sample error"}, \code{"sample accuracy"} or one of the options provided by \code{\link[ROCR]{performance}}.}
  \item{...}{Further arguments that may be used by \code{\link[ROCR]{performance}}.}
}
\details{
  Most performance metrics are provided by ROCR, so only work for two-class datasets. If
  \code{\link{runTests}} was run in resampling mode, one performance measure is produced
  for every resampling. If the leave-out mode was used, then the predictions are
  concatenated, and one performance measure is calcuated for all classifications.
  
  A variety of other metrics are also implemented in ClassifyR and are suitable for evaluating
  a multi-class classification. \code{"balanced"} calculates the balanced error rate
  and is better suited to class-imbalanced datasets than the ordinary error rate. \code{"sample error"} calculates the error rate of each sample individually. This may help to identify
  which samples are contributing the most to the error rate. \code{"sample accuracy"} causes
  the sample-wise accuracy to be computed.
}
\value{
  An updated \code{\linkS4class{ClassifyResult}} object, with new information in the \code{performance} slot.
}
\examples{
  predictTable <- data.frame(sample = 1:10,
                             label = factor(sample(LETTERS[1:2], 50, replace = TRUE)))
  actual <- factor(sample(LETTERS[1:2], 10, replace = TRUE))                             
  result <- ClassifyResult("Example", "Differential Expression", "A Selection",
                           paste("A", 1:10, sep = ''), paste("Gene", 1:50, sep = ''),
                           list(1:50, 1:50), list(1:5, 6:15),
                           list(predictTable), actual, list("leave", 2))
  result <- calcPerformance(result, "balanced") 
  performance(result)
}
\author{Dario Strbenac}
