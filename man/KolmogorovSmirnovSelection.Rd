\name{KolmogorovSmirnovSelection}
\alias{KolmogorovSmirnovSelection}
\alias{KolmogorovSmirnovSelection,matrix-method}
\alias{KolmogorovSmirnovSelection,ExpressionSet-method}
\title{Selection of Differential Distributions with Kolmogorov Smirnov Distance}
\description{Ranks features by largest Kolmogorov Smirnov distance and chooses
the features which have best resubstitution error rate.
       }
\usage{
  \S4method{KolmogorovSmirnovSelection}{matrix}(expression, classes, ...)
  \S4method{KolmogorovSmirnovSelection}{ExpressionSet}(expression, nFeatures, trainParams,
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
                 \code{\link{ExpressionSet}} method. For the \code{\link{ExpressionSet}} method,
                the options to be passed to function \code{\link{ks.test}}.}                       
       \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                      A higher number will produce more messages.}
       }
\details{
  Features are sorted in order of biggest distance to smallest. The top number of features
  is used in a classifier, to determine which number of features has the best resubstitution
  error rate.
}
\value{
  A vector of feature indices for the features that gave the best resubstitution error rate.
}
\author{Dario Strbenac}

\examples{
  if(require(sparsediscrim))
  {
    # First 25 samples are mixtures of two normals. Last 25 samples are one normal.
    genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(50, 5, 1), rnorm(50, 15, 1)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn) rnorm(100, 9, 3)))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    KolmogorovSmirnovSelection(genesMatrix, classes, nFeatures = seq(10, 100, 10),
                               trainParams = TrainParams(), predictParams = PredictParams())
  }
}
