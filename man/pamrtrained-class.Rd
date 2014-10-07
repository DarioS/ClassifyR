\name{pamrtrained}
\docType{class}
\alias{pamrtrained}
\alias{pamrtrained-class}

\title{Trained pamr Object}
\description{
  Enables dispatching on it.
}

\section{Summary}{
  A method which summarises the results is available.
  \code{result} is a \code{ClassifyResult} object.

  \describe{
    \item{}{
      \code{show(result)}{Prints a short summary of what \code{result} contains.}
    }}
}

\author{Dario Strbenac}
\examples{
  genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(100, 9, 1)))
  genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn)
                                 c(rnorm(75, 9, 1), rnorm(25, 14, 1))))
  classes <- factor(rep(c("Poor", "Good"), each = 25))
    
  result <- nearestShrunkenCentroidTrainInterface(genesMatrix, classes)
  class(result)
}
