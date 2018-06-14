\name{SelectResult}
\docType{class}
\alias{SelectResult}
\alias{SelectResult-class}
\alias{SelectResult,character,character,list,list-method}
\alias{show,SelectResult-method}

\title{Container for Storing Feature Selection Results}
\description{
  Contains a list of ranked indices or names of features, from most discriminative to least discriminative,
  and a list of indicies of features selected for use in classification. The names or indices will be in a data frame
  if the input dataset is a \code{\link{MultiAssayExperiment}}, with the first column containing the name of the
  data table the feature is from and the second column the index or name of the feature. Each vector or data frame element
  in the list corresponds to a particular iteration of classifier training. Nested lists will be present if
  the permutation and folding cross-validation scheme was used. This class is not intended to be created by the user,
  but could be used in another software package.
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
      \code{show(result)}{ Prints a short summary of what \code{result} contains.}
    }}
}

\author{Dario Strbenac}
\examples{
    SelectResult("Asthma", "Moderated t-test", list(1:50), list(1:10))
}
