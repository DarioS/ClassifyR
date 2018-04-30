\name{functionOrNULL}
\docType{class}
\alias{functionOrNULL}
\alias{functionOrNULL-class}

\title{Union of A Function and NULL}
\description{
  Allows a slot to be either a function or empty.
}

\author{Dario Strbenac}
\examples{
  PredictParams(NULL, getClasses = function(result) result)
  PredictParams(predict, getClasses = function(result) result[["classes"]])
}