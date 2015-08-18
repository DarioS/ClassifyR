\name{likelihoodRatioSelection}
\alias{likelihoodRatioSelection}
\alias{likelihoodRatioSelection,matrix-method}
\alias{likelihoodRatioSelection,ExpressionSet-method}
\title{Selection of Differential Distributions with Likelihood Ratio Statistic}
\description{Ranks features by largest ratio and chooses the features which have
             the best resubstitution performance.
       }
\usage{
  \S4method{likelihoodRatioSelection}{matrix}(expression, classes, ...)
  \S4method{likelihoodRatioSelection}{ExpressionSet}(expression, datasetName, trainParams, predictParams,
                                      resubstituteParams, alternative = c(location = "different", scale = "different"),
                                      ..., selectionName = "Likelihood Ratio Test (Normal)", verbose = 3)
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
  \item{alternative}{A vector of length 2. The first element specifies the location of
                    the alternate hypothesis. The second element specifies the scale of
                    the alternate hypothesis. Acceptable values are "same" or "different".}
  \item{...}{Either variables passed from the \code{\link{matrix}} method to the
             \code{\link{ExpressionSet}} method or variables passed to \code{\link{getLocationsAndScales}}
             from the \code{\link{ExpressionSet}} method.}
  \item{selectionName}{A name to identify this selection method by. Stored in the result.}
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  Likelihood ratio test of null hypothesis that the location and scale are the same for
  both groups, and an alternate hypothesis that is specified by parameters. The location and scale
  of features is calucated by \code{\link{getLocationsAndScales}}. The distribution fitted in the 
  normal distribution.
}
\value{
  A list of length 2. The first element has the features ranked from most important to least important.
  The second element has the features that were selected to be used for classification.
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
    likelihoodRatioSelection(genesMatrix, classes, "Example",
                             trainParams = TrainParams(naiveBayesKernel, FALSE, TRUE),
                             predictParams = PredictParams(function(){}, FALSE, getClasses = function(result) result),
                             resubstituteParams = ResubstituteParams(nFeatures = seq(10, 100, 10), performanceType = "balanced", better = "lower"))
  }
}
