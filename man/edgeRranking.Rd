\name{edgeRranking}
\alias{edgeRranking}
\alias{edgeRranking,matrix-method}
\alias{edgeRranking,DataFrame-method}
\alias{edgeRranking,MultiAssayExperiment-method}
\title{Feature Ranking Based on Differential Expression for Count Data}
\description{Performs a differential expression analysis between classes and ranks
             the features based on test statistics. The data may have overdispersion and
             this is modelled.
}
\usage{
  \S4method{edgeRranking}{matrix}(counts, classes, ...)
  \S4method{edgeRranking}{DataFrame}(counts, classes,
                   normFactorsOptions = NULL, dispOptions = NULL, fitOptions = NULL,
               verbose = 3)
  \S4method{edgeRranking}{MultiAssayExperiment}(counts, targets = NULL, ...)
}
\arguments{
  \item{counts}{Either a \code{\link{matrix}} or \code{\link{MultiAssayExperiment}} containing the
                unnormalised counts.}
  \item{classes}{A vector of class labels of class \code{\link{factor}} of the same length
                 as the number of samples in \code{measurements}. Not used if \code{measurements}
                 is a \code{MultiAssayExperiment} object.}
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables of counts to be used.}              
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method.}
  \item{normFactorsOptions}{A named \code{list} of any options to be passed to \code{\link[edgeR]{calcNormFactors}}.}
  \item{dispOptions}{A named \code{list} of any options to be passed to \code{\link[edgeR]{estimateDisp}}.}
  \item{fitOptions}{A named \code{list} of any options to be passed to \code{\link[edgeR]{glmFit}}.}                
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  The differential expression analysis follows the standard \code{\link[edgeR]{edgeR}}
  steps of estimating library size normalisation factors, calculating dispersion,
  in this case robustly, and then fitting a generalised linear model followed by
  a likelihood ratio test.
  
  Data tables which consist entirely of non-numeric data cannot be analysed. If \code{measurements}
  is an object of class \code{MultiAssayExperiment}, the factor of sample classes must be stored
  in the DataFrame accessible by the \code{colData} function with column name \code{"class"}.  
}
\value{
  A vector or data frame (if \code{MultiAssayExperiment} input) of features, from the most promising features in the first position to the least promising feature in the last position.
}
\author{Dario Strbenac}
\references{
  edgeR: a Bioconductor package for differential expression analysis of digital gene
  expression data, Mark D. Robinson, Davis McCarthy, and Gordon Smyth, 2010,
  \emph{Bioinformatics}, Volume 26 Issue 1, \url{https://academic.oup.com/bioinformatics/article/26/1/139/182458}.
}

\examples{
  if(require(parathyroidSE) && require(PoiClaClu))
  {
    data(parathyroidGenesSE)
    expression <- assays(parathyroidGenesSE)[[1]]
    sampleNames <- paste("Sample", 1:ncol(parathyroidGenesSE))
    colnames(expression) <- sampleNames
    DPN <- which(colData(parathyroidGenesSE)[, "treatment"] == "DPN")
    control <- which(colData(parathyroidGenesSE)[, "treatment"] == "Control")
    expression <- expression[, c(control, DPN)]
    classes <- factor(rep(c("Contol", "DPN"), c(length(control), length(DPN))))
    expression <- expression[rowSums(expression > 1000) > 8, ] # Make small data set.
    
    ranked <- edgeRranking(expression, classes)
                                        
    head(ranked)
    plotFeatureClasses(expression, classes, "ENSG00000044574",
                       dotBinWidth = 500, xAxisLabel = "Unnormalised Counts")
  }
}
