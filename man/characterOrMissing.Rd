\name{characterOrMissing}
\docType{class}
\alias{characterOrMissing}
\alias{characterOrMissing-class}

\title{Union of character and missing}
\description{
  Allows a method parameter to be either a character or missing.
}

\author{Dario Strbenac}
\examples{
  ModellingParams("upsample")
  ModellingParams() # Default is "downsample" if it is missing.
}
