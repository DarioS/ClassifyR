\name{ResubstituteParams}
\docType{class}
\alias{ResubstituteParams}
\alias{ResubstituteParams-class}
\alias{ResubstituteParams,missing,ANY-method}
\alias{ResubstituteParams,numeric,character-method}

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
      \code{ResubstituteParams(nFeatures, performanceType)}
      Creates a ResubstituteParams object, storing information about the number of top
      features to calculate the performance measure for and the performance measure to use.
  \describe{
            \item{\code{nFeatures}}{A vector for the top number of features to test the
                       resubstitution error rate for.}
	          \item{\code{performanceType}}{One of the eleven types of performance metrics which can be
	                                        calculated by \code{\link{calcCVperformance}}.}
            \item{\code{intermediate}}{Character vector. Names of any variables created in prior stages by
            \code{\link{runTest}} that need to be passed to the classifier.}            
            \item{\code{...}}{Other named parameters which will be used by the classifier.}
	         }
    }
  }
}
\author{Dario Strbenac}
\examples{
  ResubstituteParams(nFeatures = seq(25, 1000, 25), performanceType = "error")
}