\name{FeatureSetCollection}
\docType{class}
\alias{FeatureSetCollection}
\alias{FeatureSetCollection-class}
\alias{FeatureSetCollection,list-method}
\alias{length,FeatureSetCollection-method}
\alias{show,FeatureSetCollection-method}
\alias{[,FeatureSetCollection,integerOrNumeric,missing,ANY-method}
\alias{[[,FeatureSetCollection,ANY,missing-method}

\title{Container for Storing A Collection of Sets}
\description{
  This container is the required storage format for a collection of sets. Typically, the elements of
  a set will either be a set of proteins (i.e. character vector) which perform a particular biological process
  or a set of binary interactions (i.e. Two-column matrix of feature identifiers).
}
\section{Constructor}{
  \describe{
    \item{}{
      \code{FeatureSetCollection(sets)}}
    }
  \describe{
            \item{\code{sets}}{A named list. The names of the list describe the sets and the elements
                               of the list specify the features which comprise the sets.}
           }
}

\section{Summary}{
  A method which summarises the results is available.
  \code{featureSets} is a \code{FeatureSetCollection} object.

  \describe{
    \item{}{
      \code{show(featureSets)}{Prints a short summary of what \code{featureSets} contains.}
    }}
}

\section{Subsetting}{
  The FeatureSetCollection may be subsetted or a single set may be extracted as a vector.
  \code{featureSets} is a \code{FeatureSetCollection} object.
  
  \describe{
    \item{}{
      \code{featureSets[i:j]}:
      Reduces the object to a subset of the feature sets between elements \code{i} and \code{j}
      of the collection.
    }
    \item{}{
      \code{featureSets[[i]]}:
      Extract the feature set identified by \code{i}. \code{i} may be a numeric index
      or the character name of a feature set.
    }    
  }
}

\author{Dario Strbenac}
\examples{
    ontology <- list(c("SESN1", "PRDX1", "PRDX2", "PRDX3", "PRDX4", "PRDX5", "PRDX6",
                       "LRRK2", "PARK7"),
                     c("ATP7A", "CCS", "NQO1", "PARK7", "SOD1", "SOD2", "SOD3",
                       "SZT2", "TNF"),
                     c("AARS", "AIMP2", "CARS", "GARS", "KARS", "NARS", "NARS2",
                       "LARS2", "NARS", "NARS2", "RGN", "UBA7"),
                     c("CRY1", "CRY2", "ONP1SW", "OPN4", "RGR"),
                     c("ESRRG", "RARA", "RARB", "RARG", "RXRA", "RXRB", "RXRG"),
                     c("CD36", "CD47", "F2", "SDC4"),
                     c("BUD31", "PARK7", "RWDD1", "TAF1")
                     )
    names(ontology) <- c("Peroxiredoxin Activity", "Superoxide Dismutase Activity",
                         "Ligase Activity", "Photoreceptor Activity",
                         "Retinoic Acid Receptor Activity",
                         "Thrombospondin Receptor Activity",
                         "Regulation of Androgen Receptor Activity")
                         
    featureSets <- FeatureSetCollection(ontology)
    featureSets
    featureSets[3:5]
    featureSets[["Photoreceptor Activity"]]
    
    subNetworks <- list(MAPK = matrix(c("NRAS", "NRAS", "NRAS", "BRAF", "MEK",
                                        "ARAF", "BRAF", "CRAF", "MEK", "ERK"), ncol = 2),
                        P53 = matrix(c("ATM", "ATR", "ATR", "P53",
                                       "CHK2", "CHK1", "P53", "MDM2"), ncol = 2)
                        )
    networkSets <- FeatureSetCollection(subNetworks)                        
    networkSets
}
