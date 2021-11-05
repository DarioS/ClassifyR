\name{TransformParams}
\docType{class}
\alias{TransformParams}
\alias{TransformParams-class}
\alias{TransformParams,ANY-method}
\alias{TransformParams,function-method}
\alias{show,TransformParams-method}

\title{Parameters for Data Transformation}
\description{
  Collects and checks necessary parameters required for transformation. The 
  empty constructor is for when no data transformation is desired. See
  \code{\link{subtractFromLocation}} for an example of such a function.
}
\section{Constructor}{
  \describe{
    \item{}{
      \code{TransformParams(transform, characteristics = DataFrame(),
      intermediate = character(0), ...)}
      Creates a TransformParams object which stores the function which will do the
      transformation and parameters that the function will use.
	\describe{
            \item{\code{transform}}{A function which will do the transformation. The first
            argument must be a \code{\link{DataFrame}} object.}
            \item{\code{characteristics}}{A \code{\link{DataFrame}} describing the characteristics
            of data transformation to be done. First column must be named \code{"charateristic"} and
            second column must be named \code{"value"}. If using wrapper functions for data transformation
            in this package, the data transformation name will automatically be generated and therefore it is
            not necessary to specify it.}
            \item{\code{intermediate}}{Character vector. Names of any variables created in prior stages by
            \code{\link{runTest}} that need to be passed to a feature selection function.}            
	          \item{\code{...}}{Other named parameters which will be used by the transformation function.}
	         }
    }
  }
}

\section{Summary}{
  A method which summarises the parameter settings.
  \code{transformParams} is a \code{TransformParams} object.

  \describe{
    \item{}{
      \code{show(transformParams)}: Prints a short summary of what \code{transformParams} contains.
    }}
}

\author{Dario Strbenac}
\examples{
  transformParams <- TransformParams(subtractFromLocation, location = "median")
  # Subtract all values from training set median, to obtain absolute deviations.
}