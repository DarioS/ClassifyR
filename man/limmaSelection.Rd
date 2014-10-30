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
  \S4method{limmaSelection}{ExpressionSet}(expression, nFeatures, trainParams,
                            predictParams, ..., verbose = 3)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    the training data. For a matrix, the rows are features, and the columns
                    are samples.}
       \item{classes}{A vector of class labels.}
       \item{nFeatures}{A vector for the top number of features to test the resubstitution error
                       for.}
       \item{trainParams}{A container of class \code{\link{TrainParams}} describing the
                     classifier to use for training.}
       \item{predictParams}{A container of class \code{\link{PredictParams}} describing how
                       prediction is to be done.}
       \item{...}{For the \code{\link{matrix}} method, variables passed to the
                 \code{\link{ExpressionSet}} method. For the \code{\link{ExpressionSet}}
                 method, extra parameters that are passed to \code{\link{TrainParams}}}.
       \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                      This function only prints progress messages if the value is 3.}
       }
\details{
  This selection method looks for differential expression. It uses a moderated t-test.
}
\value{
  A list of length 2. The first element has the features ranked from most important to least important.
  The second element has the features that were selected to be used for classification.
}

\author{Dario Strbenac}
\references{
Limma: linear models for microarray data, Gordon Smyth, 2005, In: Bioinformatics and Computational Biology Solutions using R and Bioconductor, Springer, New York, pages 397-420.
}

\examples{
  if(require(sparsediscrim))
  {
    # Samples in one class with differential expression to other class.
    genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(100, 9, 1)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn)
                                 c(rnorm(75, 9, 1), rnorm(25, 14, 1))))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    
    limmaSelection(genesMatrix, classes, nFeatures = seq(10, 100, 10),
                    trainParams = TrainParams(), predictParams = PredictParams())
  }
}
