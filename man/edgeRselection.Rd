\name{edgeRselection}
\alias{edgeRselection}
\alias{edgeRselection,matrix-method}
\alias{edgeRselection,ExpressionSet-method}
\title{Feature Selection Based on Differential Expression for RNA-seq}
\description{Performs a differential expression analysis between classes and chooses
the features which have best resubstitution error rate.
}
\usage{
  \S4method{edgeRselection}{matrix}(expression, classes, ...)
  \S4method{edgeRselection}{ExpressionSet}(expression, nFeatures, normFactorsOptions = NULL,
                           dispOptions = NULL, fitOptions = NULL, trainParams,
                           predictParams, verbose = 3)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing the
                    expression values.}
  \item{classes}{A vector of class labels.}
  \item{nFeatures}{A vector for the top number of features to test the resubstitution error
                   for.}
  \item{normFactorsOptions}{A named \code{\link{list}} of any options to be passed to
                            \code{\link{calcNormFactors}}.}
  \item{dispOptions}{A named \code{\link{list}} of any options to be passed to
                            \code{\link{estimateDisp}}.}
  \item{fitOptions}{A named \code{\link{list}} of any options to be passed to
                            \code{\link{glmFit}}.}                
  \item{trainParams}{A container of class \code{\link{TrainParams}} describing the
                     classifier to use for training.}
  \item{predictParams}{A container of class \code{\link{PredictParams}} describing how
                       prediction is to be done.}
  \item{...}{Unused variables from the \code{\link{matrix}} method passed
                       to the \code{\link{ExpressionSet}} method.}
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 A higher number will produce more messages.}
}
\details{
  The differential expression analysis follows the standard \code{\link{edgeR}}
  steps of estimating library size normalisation factors, calculating dispersion,
  in this case robustly, and then fitting a generalised linear model followed by
  a likelihood ratio test.
}
\value{
  A vector of feature indices for the features that gave the best resubstitution error rate.
}
\author{Dario Strbenac}
\references{
  edgeR: a Bioconductor package for differential expression analysis of digital gene
  expression data, Mark D. Robinson, Davis McCarthy, and Gordon Smyth, 2010,
  Bioinformatics, Volume 26 Issue 1, \url{bioinformatics.oxfordjournals.org/content/26/1/139}.
}

\examples{
  if(require(parathyroidSE) && require(sparsediscrim) && require(PoiClaClu))
  {
    data(parathyroidGenesSE)
    expression <- assays(parathyroidGenesSE)[[1]]
    DPN <- which(colData(parathyroidGenesSE)[, "treatment"] == "DPN")
    control <- which(colData(parathyroidGenesSE)[, "treatment"] == "Control")
    expression <- expression[, c(control, DPN)]
    classes <- rep(c("Contol", "DPN"), c(length(control), length(DPN)))
    expression <- expression[rowSums(expression > 10) > 8, ]
    selected <- edgeRselection(expression, classes, nFeatures = seq(100, 1000, 100),
                               trainParams = TrainParams(classify, TRUE, TRUE),
                               predictParams = PredictParams(function(){}, FALSE, 
                                                             function(result) result[["ytehat"]]))
  }
}
