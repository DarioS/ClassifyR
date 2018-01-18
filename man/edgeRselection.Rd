\name{edgeRselection}
\alias{edgeRselection}
\alias{edgeRselection,matrix-method}
\alias{edgeRselection,MultiAssayExperiment-method}
\title{Feature Selection Based on Differential Expression for Count Data}
\description{Performs a differential expression analysis between classes and chooses
             the features which have best resubstitution performance. The data may have
             overdispersion and this is modelled.
}
\usage{
  \S4method{edgeRselection}{matrix}(counts, classes, ...)
  \S4method{edgeRselection}{MultiAssayExperiment}(counts, datasetName, normFactorsOptions = NULL,
                           dispOptions = NULL, fitOptions = NULL, trainParams,
                           predictParams, resubstituteParams, selectionName = "edgeR LRT", verbose = 3)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing the
                    expression values.}
  \item{classes}{A vector of class labels.}
  \item{...}{Unused variables from the \code{matrix} method passed
                       to the \code{ExpressionSet} method.}  
  \item{datasetName}{A name for the dataset used. Stored in the result.}                  
  \item{normFactorsOptions}{A named \code{list} of any options to be passed to
                            \code{\link[edgeR]{calcNormFactors}}.}
  \item{dispOptions}{A named \code{list} of any options to be passed to
                            \code{\link[edgeR]{estimateDisp}}.}
  \item{fitOptions}{A named \code{list} of any options to be passed to
                            \code{\link[edgeR]{glmFit}}.}                
  \item{trainParams}{A container of class \code{\link{TrainParams}} describing the
                     classifier to use for training.}
  \item{predictParams}{A container of class \code{\link{PredictParams}} describing how
                       prediction is to be done.}
  \item{resubstituteParams}{An object of class \code{\link{ResubstituteParams}}
                       describing the performance measure to consider and the numbers of
                       top features to try for resubstitution classification.}                   
  \item{selectionName}{A name to identify this selection method by. Stored in the result.}
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  The differential expression analysis follows the standard \code{\link[edgeR]{edgeR}}
  steps of estimating library size normalisation factors, calculating dispersion,
  in this case robustly, and then fitting a generalised linear model followed by
  a likelihood ratio test.
}
\value{
  An object of class \code{\link{SelectResult}} or a list of such objects, if the classifier which was used
  for determining resubstitution error rate made a number of prediction varieties.
}
\author{Dario Strbenac}
\references{
  edgeR: a Bioconductor package for differential expression analysis of digital gene
  expression data, Mark D. Robinson, Davis McCarthy, and Gordon Smyth, 2010,
  Bioinformatics, Volume 26 Issue 1, \url{https://academic.oup.com/bioinformatics/article/26/1/139/182458}.
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
    names(classes) <- colnames(expression)
    expression <- expression[rowSums(expression > 1000) > 8, ] # Make small dataset.
    
    getClasses <- function(result) result[["ytehat"]]
    selection <- edgeRselection(expression, classes, "DPN Treatment",
                   trainParams = TrainParams(classifyInterface, doesTests = TRUE),
                   predictParams = PredictParams(function(){}, getClasses = getClasses),
                   resubstituteParams = ResubstituteParams(nFeatures = seq(10, 100, 10),
                                        performanceType = "balanced", better = "lower"))
                                        
    head(selection@rankedFeatures)
    plotFeatureClasses(expression, classes, 77, dotBinWidth = 500,
                       xAxisLabel = "Unnormalised Counts")
  }
}
