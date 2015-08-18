\name{KolmogorovSmirnovSelection}
\alias{KolmogorovSmirnovSelection}
\alias{KolmogorovSmirnovSelection,matrix-method}
\alias{KolmogorovSmirnovSelection,ExpressionSet-method}
\title{Selection of Differential Distributions with Kolmogorov-Smirnov Distance}
\description{Ranks features by largest Kolmogorov-Smirnov distance and chooses
the features which have best resubstitution performance.
       }
\usage{
  \S4method{KolmogorovSmirnovSelection}{matrix}(expression, classes, ...)
  \S4method{KolmogorovSmirnovSelection}{ExpressionSet}(expression, datasetName, trainParams,
                                       predictParams, resubstituteParams, ...,
                                       selectionName, verbose = 3)
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
             method. For the \code{\link{ExpressionSet}} method, the options to be passed to
             function \code{\link{ks.test}}.}
  \item{selectionName}{A name to identify this selection method by. Stored in the result.}             
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
               This function only prints progress messages if the value is 3.}
}
\details{
  Features are sorted in order of biggest distance to smallest. The top number of features
  is used in a classifier, to determine which number of features has the best resubstitution
  performance.
}
\value{
  An object of class \code{\link{SelectResult}} or a list of such objects, if the classifier which was used
  for determining resubstitution error rate made a number of prediction varieties.
}
\author{Dario Strbenac}

\examples{
  if(require(sparsediscrim))
  {
    # First 20 features have bimodal distribution for Poor class. Other 80 features have normal distribution for
    # both classes.
    genesMatrix <- sapply(1:25, function(sample) c(rnorm(20, sample(c(8, 12), 20, replace = TRUE), 1), rnorm(80, 10, 1)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample) rnorm(100, 10, 1)))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    KolmogorovSmirnovSelection(genesMatrix, classes, "Example",
                               trainParams = TrainParams(naiveBayesKernel, FALSE, doesTests = TRUE),
                               predictParams = PredictParams(function(){}, FALSE, getClasses = function(result) result),
                               resubstituteParams = ResubstituteParams(nFeatures = seq(10, 100, 10), performanceType = "balanced", better = "lower"))
  }
}
