\name{listOrCharacterOrNULL}
\docType{class}
\alias{listOrCharacterOrNULL}
\alias{listOrCharacterOrNULL-class}

\title{Union of a List and a Character Vector and NULL}
\description{
  Allows a slot to be either a list, character or a NULL.
}

\author{Dario Strbenac}
\examples{
  setClass("HardClassifier", representation(model = "listOrCharacterOrNULL"))
  classifier <- new("HardClassifier", model = "Good") # Optimistic classifier.
}
