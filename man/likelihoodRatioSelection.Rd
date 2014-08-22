\name{likelihoodRatioSelection}
\alias{likelihoodRatioSelection}
\alias{likelihoodRatioSelection,matrix-method}
\alias{likelihoodRatioSelection,ExpressionSet-method}
\title{Selection of Differential Distributions with Likelihood Ratio Statistic}
\description{Ranks features by largest ratio and chooses
the features which have best resubstitution error rate.
       }
\usage{
  \S4method{likelihoodRatioSelection}{matrix}(expression, classes, ...)
  \S4method{likelihoodRatioSelection}{ExpressionSet}(expression, nFeatures, trainParams, predictParams,
                                      alternative = c(location = "different", scale = "different"),
                                      ..., verbose = 3)
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
       \item{alternative}{A vector of length 2. The first element specifies the location of
                          the alternate hypothesis. The second element specifies the scale of
                          the alternate hypothesis. Acceptable values are "same" or "different".}
       \item{...}{Variables passed to \code{\link{getLocationsAndScales}}.}                      
       \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                      A higher number will produce more messages.}
       }
\details{
  Likelihood ratio test of null hypothesis that the location and scale are the same for
  both groups, and an alternate hypothesis that is specified by parameters. The location and scale
  of features is calucated by \code{\link{getLocationsAndScales}}.
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
    likelihoodRatioSelection(genesMatrix, classes, nFeatures = seq(10, 100, 10),
                    trainParams = TrainParams(), predictParams = PredictParams())
  }
}
