\name{medianDifferenceSelection}
\alias{medianDifferenceSelection}
\alias{medianDifferenceSelection,matrix-method}
\alias{medianDifferenceSelection,ExpressionSet-method}
\title{Selection of Differential Expression by Comparing Differences in Medians of Groups}
\description{Ranks features by largest absolute difference of group medians and chooses
  the features which have best resubstitution performance.
}
\usage{
  \S4method{medianDifferenceSelection}{matrix}(expression, classes, ...)
  \S4method{medianDifferenceSelection}{ExpressionSet}(expression, datasetName,
                                            trainParams, predictParams, resubstituteParams,
                                            selectionName = "Difference of Group Medians", verbose = 3)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
    the training data. For a matrix, the rows are features, and the columns
    are samples.}
  \item{classes}{A vector of class labels.}
  \item{...}{For the \code{\link{matrix}} method, variables passed to the
    \code{\link{ExpressionSet}} method.}  
  \item{datasetName}{A name for the dataset used. Stored in the result.}    
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
\value{
  An object of class \code{\link{SelectResult}} or a list of such objects, if the classifier which was used
  for determining resubstitution error rate made a number of prediction varieties.
}
\author{Dario Strbenac}

\examples{
  if(require(sparsediscrim))
  {
    # Genes 76 to 100 have differential expression.
    genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
                                 c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    
    medianDifferenceSelection(genesMatrix, classes, datasetName = "Example",
                              trainParams = TrainParams(), predictParams = PredictParams(),
                              resubstituteParams = ResubstituteParams(nFeatures = seq(10, 100, 10),
                                                                      performanceType = "balanced", better = "lower"))
  }
}
