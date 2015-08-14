\name{functionOrList}
\docType{class}
\alias{functionOrList}
\alias{functionOrList-class}

\title{Union of Functions and List of Functions}
\description{
  Allows a slot to be either a function or a list of functions.
}

\author{Dario Strbenac}
\examples{
  SelectParams(limmaSelection)
  SelectParams(list(limmaSelection, leveneSelection), "Ensemble Selection")
}
