\name{SelectionParams}
\docType{class}
\alias{SelectionParams}
\alias{SelectionParams-class}
\alias{SelectionParams,ANY-method}
\alias{SelectionParams,functionOrList-method}

\title{Parameters for Feature Selection}
\description{
  Collects and checks necessary parameters required for feature selection. The 
  empty constructor is provided for convenience.
}
\section{Constructor}{
  \describe{
    \item{}{
      \code{SelectionParams()}
      Creates a default SelectionParams object. This uses a limma t-test and tries 100,
      200, 300, 400, 500 features, and picks the number of features with the best
      resubstitution error rate. Users should create an appropriate \code{SelectionParams}
      object for the characteristics of their data, once they are familiar with this software.
  }
    \item{}{
      \code{SelectionParams(featureSelection, minPresence = 1, intermediate = character(0),
                   subsetExpressionData = TRUE, ...)}
      Creates a SelectionParams object which stores the function which will do the
      selection and parameters that the function will use.
	\describe{
            \item{\code{featureSelection}}{Either a function which will do the selection
            or a list of such functions. For a particular function, the first argument must
            be an \code{\link[Biobase]{ExpressionSet}} object. The function's
            return value must be a vector of row indices of genes that were selected.}
            \item{\code{minPresence}}{If a list of functions was provided, how many of those
            must a feature have been selected by to be used in classification. 1 is equivalent
            to a set union and a number the same length as \code{featureSelection}
            is equivalent to set intersection.}
            \item{\code{intermediate}}{Character vector. Names of any variables created in prior stages by
            \code{\link{runTest}} that need to be passed to a feature selection function.}
            \item{\code{subsetExpressionData}}{Whether to subset the expression data, after
            selection has been done.}
	          \item{\code{...}}{Other named parameters which will be used by the selection function.
            If \code{featureSelection} was a list of functions, this must be a list of lists, as
            long as \code{featureSelection}.}
	         }
    }
  }
}

\author{Dario Strbenac}
\examples{
  SelectionParams(limmaSelection, nFeatures = c(25, 50, 75, seq(100, 1000, 100)))
  
  # For pamr shrinkage selection.
  SelectionParams(nearestShrunkenCentroidSelectionInterface, intermediate = "trained",
                  subsetExpressionData = FALSE)
}