\name{leveneSelection}
\alias{leveneSelection}
\alias{leveneSelection,matrix-method}
\alias{leveneSelection,ExpressionSet-method}
\title{Selection of Differential Variability with Levene Statistic}
\description{Ranks features by largest Levene statistic and chooses
the features which have best resubstitution performance.
       }
\usage{
  \S4method{leveneSelection}{matrix}(expression, classes, ...)
  \S4method{leveneSelection}{ExpressionSet}(expression,
                             trainParams, predictParams, resubstituteParams, verbose = 3)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    the training data. For a matrix, the rows are features, and the columns
                    are samples.}
  \item{classes}{A vector of class labels.}
  \item{trainParams}{A container of class \code{\link{TrainParams}} describing the
                     classifier to use for training.}
  \item{predictParams}{A container of class \code{\link{PredictParams}} describing how
                       prediction is to be done.}
  \item{resubstituteParams}{An object of class \code{\link{ResubstituteParams}}
                            describing the performance measure to consider and the numbers of
                            top features to try for resubstitution classification.}                               
  \item{...}{For the \code{\link{matrix}} method, variables passed to the
                 \code{\link{ExpressionSet}} method.}                     
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  Levene's statistic for unequal variance between groups is a robust version of
  Bartlett's statistic.
}
\value{
  A list of length 2. The first element has the features ranked from most important to least important.
  The second element has the features that were selected to be used for classification.
}
\author{Dario Strbenac}

\examples{
  if(require(sparsediscrim))
  {
    # Samples in one class with differential variability to other class.
    genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(100, 9, 1)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn) rnorm(100, 9, 4)))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    leveneSelection(genesMatrix, classes,
                    trainParams = TrainParams(), predictParams = PredictParams(),
                    resubstituteParams = ResubstituteParams(nFeatures = seq(10, 100, 10),
                                         performanceType = "balanced", better = "lower"))
  }
}
