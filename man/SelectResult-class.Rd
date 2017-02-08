\name{SelectResult}
\docType{class}
\alias{SelectResult}
\alias{SelectResult-class}
\alias{SelectResult,character,character,list,list-method}
\alias{show,SelectResult-method}

\title{Container for Storing Feature Selection Results}
\description{
  Contains the ranked indices or names of features, from most discriminative to least discriminative and a list
  of indicies of feature selected for use in classification. This class is not intended to be created by
  the user, but could be used in another package.
}
\section{Constructor}{
  \describe{
    \item{}{
      \code{SelectResult(datasetName, selectionName, rankedFeatures, chosenFeatures)}}
    }
  \describe{
            \item{\code{datasetName}}{A name associated with the dataset used.}
            \item{\code{selectionName}}{A name associated with the classification.}
            \item{\code{rankedFeatures}}{Indices or names of all features, from most to least
                                         discriminative.}
            \item{\code{chosenFeatures}}{Indices or names of features selected at each fold.}
           }
}

\section{Summary}{
  A method which summarises the results is available.
  \code{result} is a \code{SelectResult} object.

  \describe{
    \item{}{
      \code{show(result)}{Prints a short summary of what \code{result} contains.}
    }}
}

\author{Dario Strbenac}
\examples{
    SelectResult("Melanoma", "Moderated t-test", list(1:50), list(1:10))
}
