\name{TransformParams}
\docType{class}
\alias{TransformParams}
\alias{TransformParams-class}
\alias{TransformParams,ANY-method}
\alias{TransformParams,function-method}

\title{Parameters for Data Transformation}
\description{
  Collects and checks necessary parameters required for transformation. The 
  empty constructor is for when no data transformation is desired. One data
  transformation function is distributed. See \code{\link{subtractFromLocation}}.
}
\section{Constructor}{
  \describe{
    \item{}{
      \code{TransformParams(transform, intermediate = character(0), ...)}
      Creates a TransformParams object which stores the function which will do the
      transformation and parameters that the function will use.
	\describe{
            \item{\code{transform}}{A function which will do the transformation. The first
            argument must be an \code{\link[Biobase]{ExpressionSet}} object.}
            \item{\code{intermediate}}{Character vector. Names of any variables created in prior stages by
            \code{\link{runTest}} that need to be passed to a feature selection function.}            
	          \item{\code{...}}{Other named parameters which will be used by the transformation function.}
	         }
    }
  }
}
\author{Dario Strbenac}
\examples{
transforParams <- TransformParams(subtractFromLocation, location = "median")
# Subtract all values from training set median, to obtain absolute deviations.
}