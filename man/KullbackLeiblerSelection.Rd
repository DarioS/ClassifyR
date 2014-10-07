\name{KullbackLeiblerSelection}
\alias{KullbackLeiblerSelection}
\alias{KullbackLeiblerSelection,matrix-method}
\alias{KullbackLeiblerSelection,ExpressionSet-method}
\title{Selection of Differential Distributions with Kullback Leibler Distance}
\description{Ranks features by largest Kullback Leibler distance and chooses
the features which have best resubstitution error rate.
       }
\usage{
  \S4method{KullbackLeiblerSelection}{matrix}(expression, classes, ...)
  \S4method{KullbackLeiblerSelection}{ExpressionSet}(expression, nFeatures, trainParams,
                                      predictParams, type = c("symmetric", "dynamic"),
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
       \item{type}{A \code{character} of length 1. Either "symmetric" or "dynamic".
                   See details for definitions.} 
       \item{...}{Variables passed to \code{\link{getLocationsAndScales}}.}                       
       \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                      This function only prints progress messages if the value is 3.}
       }
\details{
  Symmetric distance is defined as \eqn{1/2 * (location_1 - location_2)^2}{0.5 * (location1 - location2)^2 / scale1^2 + (location1 - location2)^2 / scale2^2 + scale1^2 / scale2^2 + scale2^2 / scale1^2 }
  
  Dynamic distance is defined as \eqn{|location_1 - location_2| + |scale_1 -  scale_2|}{|location1 - location2| + |scale1 -  scale2|}
  
  The subscripts denote the group which the parameter is calculated for.
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
    KullbackLeiblerSelection(genesMatrix, classes, nFeatures = seq(10, 100, 10),
                             trainParams = TrainParams(), predictParams = PredictParams())
  }
}
