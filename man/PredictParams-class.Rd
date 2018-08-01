\name{PredictParams}
\docType{class}
\alias{PredictParams}
\alias{PredictParams-class}
\alias{PredictParams,ANY-method}
\alias{PredictParams,functionOrNULL-method}

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
      \code{PredictParams(predictor, intermediate = character(0), getClasses, ...)}
      Creates a PredictParams object which stores the function which will do the
      class prediction, if required, and parameters that the function will use. If the training
      function also makes predictions, this must be set to \code{NULL}.
  \describe{
            \item{\code{predictor}}{Either \code{NULL} or a function to make predictions with. If it is
            a function, then the first argument must accept the classifier made in the training step.
            The second argument must accept a \code{\link{DataFrame}} of new data.}
            \item{\code{intermediate}}{Character vector. Names of any variables created in prior stages in
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
predictParams <- PredictParams(predictor = DLDApredictInterface,
                               getClasses = function(result) result)
# For prediction by trained object created by DLDA function.
PredictParams(predictor = NULL, getClasses = function(result) result)
# For when the training function also does prediction and directly returns the
# predictions.
}