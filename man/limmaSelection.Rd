\name{limmaSelection}
\alias{limmaSelection}
\alias{limmaSelection,matrix-method}
\alias{limmaSelection,ExpressionSet-method}
\title{Selection of Differentially Expressed Features}
\description{Uses a moderated t-test with empirical Bayes shrinkage to select
             differentially expressed features.
       }
\usage{
  \S4method{limmaSelection}{matrix}(expression, classes, ...)
  \S4method{limmaSelection}{ExpressionSet}(expression, datasetName, trainParams, predictParams,
                                           resubstituteParams, ..., selectionName = "Moderated t-test", verbose = 3)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    the training data. For a matrix, the rows are features, and the columns
                    are samples.}
  \item{classes}{A vector of class labels.}                    
  \item{datasetName}{A name for the dataset used. Stored in the result.}                    
  \item{trainParams}{A container of class \code{\link{TrainParams}} describing the
                     classifier to use for training.}
  \item{predictParams}{A container of class \code{\link{PredictParams}} describing how
                       prediction is to be done.}
  \item{resubstituteParams}{An object of class \code{\link{ResubstituteParams}}
                            describing the performance measure to consider and the numbers of
                            top features to try for resubstitution classification.}                       
  \item{...}{For the \code{\link{matrix}} method, variables passed to the \code{\link{ExpressionSet}}
             method. For the \code{\link{ExpressionSet}} method, extra parameters that are
             passed to \code{\link{[limma]lmFit}}}.
  \item{selectionName}{A name to identify this selection method by. Stored in the result.}
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  This selection method looks for differential expression. It uses a moderated t-test.
}
\value{
  An object of class \code{\link{SelectResult}} or a list of such objects, if the classifier which was used
  for determining resubstitution error rate made a number of prediction varieties.
}

\author{Dario Strbenac}
\references{
Limma: linear models for microarray data, Gordon Smyth, 2005, In: Bioinformatics and Computational Biology Solutions using R and Bioconductor, Springer, New York, pages 397-420.
}

\examples{
  if(require(sparsediscrim))
  {
    # Genes 76 to 100 have differential expression.
    genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
                                 c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    
    limmaSelection(genesMatrix, classes, "Example",
                    trainParams = TrainParams(), predictParams = PredictParams(),
                    resubstituteParams = ResubstituteParams(nFeatures = seq(10, 100, 10), performanceType = "balanced", better = "lower"))
  }
}
