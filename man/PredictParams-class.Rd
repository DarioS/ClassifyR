\name{PredictParams}
\docType{class}
\alias{PredictParams}
\alias{PredictParams-class}
\alias{PredictParams,missing-method}
\alias{PredictParams,functionOrNULL-method}
\alias{show,PredictParams-method}

\title{Parameters for Classifier Prediction}
\description{
  Collects the function to be used for making predictions and any associated parameters.
  
  The function specified must return either a factor vector of class predictions, or a numeric vector
  of scores for the second class, according to the levels of the class vector of the input data set,
  or a data frame which has two columns named class and score.
}
\section{Constructor}{
  \describe{
    \item{}{
      \code{PredictParams()}
      Creates a default PredictParams object. This assumes that the object returned by
      the classifier has a list element named \code{"class"}.
  }
    \item{}{
      \code{PredictParams(predictor, characteristics = DataFrame(), intermediate = character(0), ...)}
      Creates a PredictParams object which stores the function which will do the
      class prediction, if required, and parameters that the function will use. If the training
      function also makes predictions, this must be set to \code{NULL}.
      }
  \describe{
            \item{\code{predictor}}{Either \code{NULL} or a function to make predictions with. If it is
            a function, then the first argument must accept the classifier made in the training step.
            The second argument must accept a \code{\link{DataFrame}} of new data.}
            \item{\code{characteristics}}{A \code{\link{DataFrame}} describing the characteristics
            of the predictor function used. First column must be named \code{"charateristic"} and
            second column must be named \code{"value"}.}
            \item{\code{intermediate}}{Character vector. Names of any variables created in prior stages in
            \code{\link{runTest}} that need to be passed to the prediction function.}
            \item{\code{...}}{Other arguments that \code{predictor} may use.}
	         }
    }
  }
  
\section{Summary}{
  A method which summarises the parameter settings.
  \code{predictParams} is a \code{PredictParams} object.

  \describe{
    \item{}{
      \code{show(predictParams)}: Prints a short summary of what \code{predictParams} contains.
    }}
}

\author{Dario Strbenac}
\examples{
predictParams <- PredictParams(predictor = DLDApredictInterface)
# For prediction by trained object created by DLDA training function.
PredictParams(predictor = NULL)
# For when the training function also does prediction and directly returns the
# predictions.
}