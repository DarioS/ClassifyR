\name{numericOrMissing}
\docType{class}
\alias{numericOrMissing}
\alias{numericOrMissing-class}

\title{Union of numeric and missing}
\description{
  Allows a method parameter to be either a numeric value or missing.
}

\author{Dario Strbenac}
\examples{
  CrossValParams(permutations = 10) # numeric
  CrossValParams() # missing
}
