\name{leveneSelection}
\alias{leveneSelection}
\alias{leveneSelection,matrix-method}
\alias{leveneSelection,ExpressionSet-method}
\title{Selection of Differential Variability with Levene Statistic}
\description{Ranks features by largest Levene statistic and chooses
the features which have best resubstitution error rate.
       }
\usage{
  \S4method{leveneSelection}{matrix}(expression, classes, ...)
  \S4method{leveneSelection}{ExpressionSet}(expression, nFeatures,
                             trainParams, predictParams, verbose = 3)
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
                 \code{\link{ExpressionSet}} method.}                     
       \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                      A higher number will produce more messages.}
       }
\details{
  Levene's statistic for unequal variance between groups is a robust version of
  Bartlett's statistic.
}
\value{
  A vector of feature indices for the features that gave the best resubstitution error rate.
}
\author{Dario Strbenac}

\examples{
  if(require(sparsediscrim))
  {
    # Samples in one class with differential variability to other class.
    genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(100, 9, 1)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn) rnorm(100, 9, 4)))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    leveneSelection(genesMatrix, classes, nFeatures = seq(10, 100, 10),
                    trainParams = TrainParams(), predictParams = PredictParams())
  }
}
