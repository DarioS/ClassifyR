\name{MixModelsListsSet}
\docType{class}
\alias{MixModelsListsSet}
\alias{MixModelsListsSet-class}
\alias{MixModelsListsSet,list-method}

\title{Container for a List of Lists Containing Mixture Models}
\description{
  Stores a list of lists of trained mixture models, to prevent them being unintentionally being unlisted during cross-validation. Not intended for end-user.
}
\section{Constructor}{
  \describe{
    \item{}{
      \preformatted{MixModelsListsSet(set)}
            Creates a MixModelsListsSet object which stores the mixture models.
  \describe{
            \item{\code{set}}{A list as long as the number of classes in the data set. Each element is a list, which each element of is a mixture model trained on one feature.}
	         }
    }
  }
}
\author{Dario Strbenac}
\examples{
  if(require(Rmixmod))
  {
    mixModels <- list(Good = list(mixmodCluster(rnorm(20), nbCluster = 1:2)),
                      Poor = list(mixmodCluster(rnorm(20), nbCluster = 1:2)))
    MixModelsListsSet(mixModels)
  }
}