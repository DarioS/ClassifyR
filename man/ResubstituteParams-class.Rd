\name{ResubstituteParams}
\docType{class}
\alias{ResubstituteParams}
\alias{ResubstituteParams-class}
\alias{ResubstituteParams,ANY,ANY,ANY-method}
\alias{ResubstituteParams,numeric,character,character-method}

\title{Parameters for Resubstitution Error Calculation}
\description{
  Some feature selection functions provided in the framework use resubstitution error
  rate to choose the best number of features for classification. This class stores
  parameters related to that process.
}
\section{Constructor}{
  \describe{
    \item{}{
      \code{ResubstituteParams()}
      Creates a default ResubstituteParams object. The number of features tried is
      10, 20, 30, 40, 50, 60, 70, 80, 90, 100. The performance measure used is the
      balanced error rate.
  }
    \item{}{
      \code{ResubstituteParams(nFeatures, performanceType, better = c("lower", "higher"))}
      Creates a ResubstituteParams object, storing information about the number of top
      features to calculate the performance measure for, the performance measure to use,
      and if higher or lower values of the measure are better.
  \describe{
            \item{\code{nFeatures}}{A vector for the top number of features to test the
                       resubstitution error for.}
	          \item{\code{performanceType}}{One of the eleven types of performance metrics which can be
	                                        calculated by \code{\link{calcCVperformance}}.}
            \item{\code{better}}{Either \code{"lower"} or \code{"higher"}. Determines whether higher
                       or lower values of the performance measure are desirable.}
            \item{\code{intermediate}}{Character vector. Names of any variables created in prior stages by
            \code{\link{runTest}} that need to be passed to the classifier.}            
            \item{\code{...}}{Other named parameters which will be used by the classifier.}
	         }
    }
  }
}
\author{Dario Strbenac}
\examples{
  ResubstituteParams(nFeatures = seq(25, 1000, 25), performanceType = "error", better = "lower")
}