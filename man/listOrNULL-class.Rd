\name{listOrNULL}
\docType{class}
\alias{listOrNULL}
\alias{listOrNULL-class}

\title{Union of a List and NULL}
\description{
  Allows a slot to be either a list or a NULL.
}

\author{Dario Strbenac}
\examples{
  setClass("EasyClassifier", representation(model = "listOrNULL"))
  classifier <- new("EasyClassifier", model = NULL) # Optimistic classifier.
}
