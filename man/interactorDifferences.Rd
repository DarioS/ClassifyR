\name{interactorDifferences}
\alias{interactorDifferences}
\alias{interactorDifferences,matrix-method}
\alias{interactorDifferences,DataFrame-method}
\alias{interactorDifferences,MultiAssayExperiment-method}
\title{Convert Individual Features into Differences Between Binary Interactors Based on Known Sub-networks}
\description{This conversion is useful for creating a meta-feature table for classifier training and prediction
             based on sub-networks that were selected based on their differential correlation between classes.
}
\usage{
  \S4method{interactorDifferences}{matrix}(measurements, ...)
  \S4method{interactorDifferences}{DataFrame}(measurements, networkSets = NULL, verbose = 3)
  \S4method{interactorDifferences}{MultiAssayExperiment}(measurements, target = NULL, ...)
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
                      \code{\link{MultiAssayExperiment}} containing the training data.
                      For a \code{matrix}, the rows are features, and the columns are samples.}
  \item{networkSets}{A object of type \code{\link{FeatureSetCollection}}. The \code{sets} slot must contain
                     a list of two-column matrices with each row corresponding to a binary interaction. Such
                     sub-networks may be determined by a community detection algorithm.}
  \item{target}{If \code{measurements} is a \code{MultiAssayExperiment}, the name of the
                 data table to be used.}
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method.}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
    This function only prints progress messages if the value is 3.}
}
\details{
  The pairs of features known to interact with each other are specified by \code{networkSets}.
}
\value{
  An object of class \code{\link{DataFrame}} with one column for each interactor pair difference and one row
  for each sample. Additionally, \code{mcols(resultTable)} prodvides a \code{\link{DataFrame}} with a column
  named "original" containing the name of the sub-network each meta-feature belongs to.
}
\author{Dario Strbenac}
\references{
  Dynamic modularity in protein interaction networks predicts breast cancer outcome,
  Ian W Taylor, Rune Linding, David Warde-Farley, Yongmei Liu, Catia Pesquita, Daniel Faria,
  Shelley Bull, Tony Pawson, Quaid Morris and Jeffrey L Wrana, 2009,
  \emph{Nature Biotechnology}, Volume 27 Issue 2, \url{https://www.nature.com/articles/nbt.1522}.
}
\examples{
  networksList <- list(`A Hub` = matrix(c('A', 'A', 'A', 'B', 'C', 'D'), ncol = 2),
                       `G Hub` = matrix(c('G', 'G', 'G', 'H', 'I', 'J'), ncol = 2))
  netSets <- FeatureSetCollection(networksList)
                           
  # Differential correlation for sub-network with hub A.                                           
  measurements <- matrix(c(5.7, 10.1, 6.9, 7.7, 8.8, 9.1, 11.2, 6.4, 7.0, 5.5,
                           5.6, 9.6, 7.0, 8.4, 10.8, 12.2, 8.1, 5.7, 5.4, 12.1,
                           4.5, 9.0, 6.9, 7.0, 7.3, 6.9, 7.8, 7.9, 5.7, 8.7,
                           8.1, 10.6, 7.4, 7.15, 10.4, 6.1, 7.3, 2.7, 11.0, 9.1,
                           round(rnorm(60, 8, 1), 1)), ncol = 10, byrow = TRUE)
                         
  rownames(measurements) <- LETTERS[1:10]
  colnames(measurements) <- paste("Patient", 1:10)
  
  interactorDifferences(measurements, netSets)
}
