\name{integerOrNumeric}
\docType{class}
\alias{integerOrNumeric}
\alias{integerOrNumeric-class}

\title{Union of a Integer and a Numeric}
\description{
  Allows the same S4 subsetting function to be specified for \code{object[i]} and
  \code{object[i:j]}, where \code{i} and \code{j} are integers.
}

\author{Dario Strbenac}
\examples{
  setClass("Container", representation(scores = "numeric"))
  setMethod("[", c("Container", "integerOrNumeric", "missing", "ANY"),
           function(x, i, j, ..., drop = TRUE)
           {
             new("Container", scores = x@scores[i])
           })
           
  dataset <- new("Container", scores = 1:10)
  dataset[1] # 1 is numeric.
  dataset[4:6] # 4:6 is a sequence of integers.
}
