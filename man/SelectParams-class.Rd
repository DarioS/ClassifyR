\name{SelectParams}
\docType{class}
\alias{SelectParams}
\alias{SelectParams-class}
\alias{SelectParams,ANY-method}
\alias{SelectParams,functionOrList-method}

\title{Parameters for Feature Selection}
\description{
  Collects and checks necessary parameters required for feature selection. The 
  empty constructor is provided for convenience.
}
\section{Constructor}{
  \describe{
    \item{}{
      \code{SelectParams()}
      Creates a default SelectParams object. This uses a limma t-test and tries the top 10 to
      top 100 features in increments of 10, and picks the number of features with the best
      resubstitution balanced error rate. Users should create an appropriate \code{SelectParams}
      object for the characteristics of their data, once they are familiar with this software.
  }
    \item{}{
      \code{SelectParams(featureSelection, selectionName, minPresence = 1, intermediate = character(0),\cr
             subsetToSelections = TRUE, ...)}
      Creates a SelectParams object which stores the function which will do the
      selection and parameters that the function will use.
	\describe{
            \item{\code{featureSelection}}{Either a function which will do the selection
            or a list of such functions. For a particular function, the first argument must
            be an \code{\link{DataFrame}} object. The function's return value must be
            a \code{\link{SelectResult}} object.}
            \item{\code{selectionName}}{A name to identify this selection method by.}
            \item{\code{minPresence}}{If a list of functions was provided, how many of those
            must a feature have been selected by to be used in classification. 1 is equivalent
            to a set union and a number the same length as \code{featureSelection}
            is equivalent to set intersection.}
            \item{\code{intermediate}}{Character vector. Names of any variables created in prior stages by
            \code{\link{runTest}} that need to be passed to a feature selection function.}
            \item{\code{subsetToSelections}}{Whether to subset the data table(s), after
                                             feature selection has been done.}
	          \item{\code{...}}{Other named parameters which will be used by the selection function.
            If \code{featureSelection} was a list of functions, this must be a list of lists, as
            long as \code{featureSelection}.}
	         }
    }
  }
}

\author{Dario Strbenac}
\examples{
  if(require(sparsediscrim))
  {
    SelectParams(limmaSelection, "t-test",
                 trainParams = TrainParams(), predictParams = PredictParams(),
                 resubstituteParams = ResubstituteParams())
    
    # For pamr shrinkage selection.
    SelectParams(NSCselectionInterface, datasetName = "Cancer",
                 intermediate = "trained", subsetToSelections = FALSE)
  }
}