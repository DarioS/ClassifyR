\name{edgesToHubNetworks}
\alias{edgesToHubNetworks}

\title{
  Convert a Two-column Matrix or Data Frame into a Hub Node List
}
\description{
 Interactions between pairs of features (typically a protein-protein interaction, commonly abbreviated as PPI, database) are restructured into a named list. The name of the each element of the list is a feature and the element contains all features which have an interaction with it.
}
\usage{
  edgesToHubNetworks(edges, minCardinality = 5)
}
\arguments{
  \item{edges}{A two-column \code{matrix} or \code{data.frame} for which each row specifies a known interaction
               betwen two interactors. If feature X appears in the first column and feature Y appears in the second,
               there is no need for feature Y to appear in the first column and feature X in the second.}
  \item{minCardinality}{An integer specifying the minimum number of features to be associated with a hub feature
                        for it to be present in the result.}
}
\value{
  An object of type \code{\link{FeatureSetCollection}}.
}
\references{
  VAN: an R package for identifying biologically perturbed networks via differential variability analysis,
  Vivek Jayaswal, Sarah-Jane Schramm, Graham J Mann, Marc R Wilkins and Yee Hwa Yang, 2010,
  \emph{BMC Research Notes}, Volume 6 Article 430, \url{https://bmcresnotes.biomedcentral.com/articles/10.1186/1756-0500-6-430}.
}
\author{
  Dario Strbenac
}
\examples{
  interactor <- sapply(1:10000, function(index)
                       paste(c(sample(LETTERS, 3), sample(1:9, 1)), collapse = ''))
  otherInteractor <- sapply(1:10000, function(index)
                       paste(c(sample(LETTERS, 3), sample(1:9, 1)), collapse = ''))
  edges <- data.frame(interactor, otherInteractor, stringsAsFactors = FALSE)
  
  edgesToHubNetworks(edges, minCardinality = 4)
}

