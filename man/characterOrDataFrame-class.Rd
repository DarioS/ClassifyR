\name{characterOrDataFrame}
\docType{class}
\alias{characterOrDataFrame}
\alias{characterOrDataFrame-class}

\title{Union of a Character and a DataFrame}
\description{
  Allows a slot to be either a character or a DataFrame.
}

\author{Dario Strbenac}
\examples{
  setClass("Selections", representation(features = "characterOrDataFrame"))
  selections <- new("Selections", features = c("BRAF", "NRAS"))
  featuresTable <- DataFrame(assay = c("RNA-seq", "Mass spectrometry"),
                             feature = c("PD-1", "MITF"))
  omicsSelections <- new("Selections", features = featuresTable)
}
