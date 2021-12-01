\name{SelectParamsOrNULL}
\docType{class}
\alias{SelectParamsOrNULL}
\alias{SelectParamsOrNULL-class}

\title{Union of A SelectParams Object and NULL}
\description{
  Allows a slot to be either a SelectParams class object or empty. No constructor.
}

\author{Dario Strbenac}

\examples{
  ModellingParams(selectParams = NULL)
  ModellingParams(selectParams = SelectParams(differentMeansRanking))
}