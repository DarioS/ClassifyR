#' Asthma RNA Abundance and Patient Classes
#' 
#' Data set consists of a matrix of abundances of 2000 most variable gene
#' expression measurements for 190 samples and a factor vector of classes for
#' those samples.
#' 
#' 
#' @name asthma
#' @aliases measurements classes
#' @docType data
#' @format \code{measurements} has a row for each sample and a column for each
#' gene. \code{classes} is a factor vector with values No and Yes, indicating if
#' a partiular person has asthma or not.
#' @source A Nasal Brush-based Classifier of Asthma Identified by Machine
#' Learning Analysis of Nasal RNA Sequence Data, \emph{Scientific Reports},
#' 2018.  Webpage: \url{http://www.nature.com/articles/s41598-018-27189-4}
#' @keywords datasets
NULL


#' Human Reference Interactome
#' 
#' A collection of 45783 pairs of protein gene symbols, as determined by the
#' The Human Reference Protein Interactome Mapping Project. Self-interactions
#' have been removed.
#' 
#' 
#' @name HuRI
#' @aliases interactors
#' @docType data
#' @format \code{interactors} is a \code{\link{Pairs}} object containing each
#' pair of interacting proteins.
#' @source A Reference Map of the Human Binary Protein Interactome,
#' \emph{Nature}, 2020.  Webpage:
#' \url{http://www.interactome-atlas.org/download}
#' @keywords datasets
NULL


