\name{PredictParams}
\docType{class}
\alias{PredictParams}
\alias{PredictParams-class}
\alias{PredictParams,ANY-method}
\alias{PredictParams,function-method}

\title{Parameters for Classifier Prediction}
\description{
  Collects the function to be used for making predictions and any associated parameters.
}
\section{Constructor}{
  \describe{
    \item{}{
      \code{PredictParams()}
      Creates a default PredictParams object. This assumes that the object returned by
      the classifier has a list element named "class".
  }
    \item{}{
      \code{PredictParams(predictor, transposeExpression, multipleResults, intermediate = character(0), getClasses, ...)}
      Creates a PredictParams object which stores the function which will do the
      class prediction and parameters that the function will use.
  \describe{
            \item{\code{predictor}}{A \code{\link{function}} to make predictions with.
            The first argument must accept the classifier made in the training step.
            The second argument must accept a \code{\link{matrix}} of new data.}
            \item{\code{transposeExpression}}{Set to \code{TRUE} if \code{classifier} expects
            features as columns.}
            \item{\code{multipleResults}}{If \code{TRUE}, \code{predictor} will return a
            \code{\link{list}} of results, perhaps for different values of some parameter.
            If FALSE, a single result is returned.}
            \item{\code{intermediate}}{Character vector. Names of any variables created in prior stages by
            \code{\link{runTest}} that need to be passed to the prediction function.}               
            \item{\code{getClasses}}{A \code{\link{function}} to extract the vector of class predictions
            from the result object created by \code{predictor}.}
            \item{\code{...}}{Other arguments that \code{predictor} may use.}
	         }
    }
  }
}
\author{Dario Strbenac}

\examples{
predictParams <- PredictParams(predictor = predict, TRUE, FALSE, getClasses = function(result) result)
# For prediction by trained object created by dlda function.
PredictParams(predictor = function(){}, TRUE, FALSE, getClasses = function(result) result)
# For when the training function also does prediction and directly returns vector of predictions.
}