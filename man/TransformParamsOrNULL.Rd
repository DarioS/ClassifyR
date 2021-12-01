\name{TransformParamsOrNULL}
\docType{class}
\alias{TransformParamsOrNULL}
\alias{TransformParamsOrNULL-class}

\title{Union of A TransformParams Object and NULL}
\description{
  Allows a slot to be either a TransformParams class object or empty. No constructor.
}

\author{Dario Strbenac}

\examples{
  ModellingParams(transformParams = NULL)
  ModellingParams(transformParams = TransformParams(subtractFromLocation),
                  selectParams = SelectParams(leveneRanking))
}