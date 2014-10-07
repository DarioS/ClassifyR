\name{TrainParams}
\docType{class}
\alias{TrainParams}
\alias{TrainParams-class}
\alias{TrainParams,ANY-method}
\alias{TrainParams,function-method}

\title{Parameters for Classifier Training}
\description{
  Collects and checks necessary parameters required for classifier training. The 
  empty constructor is provided for convenience.
}
\section{Constructor}{
  \describe{
    \item{}{
      \code{TrainParams()}
      Creates a default TrainParams object. The classifier function is DLDA. Users
      should create an appropriate \code{TrainParams} object for the
      characteristics of their data, once they are familiar with this software.
  }
    \item{}{
      \code{TrainParams(classifier, transposeExpression, doesTests, ...)}
      Creates a TrainParams object which stores the function which will do the
      classifier building and parameters that the function will use.
  \describe{
            \item{\code{classifier}}{A function which will construct a classifier, and also
            possibly make the predictions. The first argument must be a \code{\link{matrix}}
            object. The second argument must be a vector of classes. The third argument must be
            \code{verbose}. If doesTests is \code{TRUE}, the third argument must be a \code{matrix}
            of test data and the fourth argument is \code{verbose}. The function's return value can
            be either a trained classifier when \code{doesTests} is \code{FALSE} or a
            vector of class predictions if \code{doesTests} is \code{TRUE}.}
	          \item{\code{transposeExpression}}{Set to \code{TRUE} if \code{classifier} expects
            features as columns.}
            \item{\code{doesTests}}{Set to \code{TRUE} if \code{classifier} also performs
            and returns predictions.}
            \item{\code{intermediate}}{Character vector. Names of any variables created in prior stages by
            \code{\link{runTest}} that need to be passed to \code{classifier}.}            
            \item{\code{...}}{Other named parameters which will be used by the classifier.}
	         }
    }
  }
}
\author{Dario Strbenac}
\examples{
if(require(sparsediscrim))
  trainParams <- TrainParams(dlda, transposeExpression = TRUE, doesTests = FALSE)
# sparsediscrim has a separate predict method for trained DLDA objects.
# dlda expects features in columns, and samples in rows.
}