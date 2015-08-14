\name{bartlettSelection}
\alias{bartlettSelection}
\alias{bartlettSelection,matrix-method}
\alias{bartlettSelection,ExpressionSet-method}
\title{Selection of Differential Variability with Bartlett Statistic}
\description{Ranks features by largest Bartlett statistic and chooses
  the features which have best resubstitution performance.
}
\usage{
  \S4method{bartlettSelection}{matrix}(expression, classes, ...)
  \S4method{bartlettSelection}{ExpressionSet}(expression, datasetName,
                                            trainParams, predictParams, resubstituteParams,
                                            selectionName = "Bartlett Test", verbose = 3)
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
\details{
  The calculation of the test statistic is performed by the \code{\link{bartlett.test}}
  function from the \code{\link{stats}} package.
}
\value{
  An object of class \code{\link{SelectResult}} or a list of such objects, if the classifier which was used
  for determining resubstitution error rate made a number of prediction varieties.
}
\author{Dario Strbenac}

\examples{
  if(require(sparsediscrim))
  {
    # Samples in one class with differential variability to other class.
    # First 20 genes are DV.
    genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 1)))
    genesMatrix <- cbind(genesMatrix, rbind(sapply(1:25, function(sample) rnorm(20, 9, 5)),
                                            sapply(1:25, function(sample) rnorm(80, 9, 1))))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    genesMatrix <- exprs(subtractFromLocation(genesMatrix, 1:ncol(genesMatrix)))
    bartlettSelection(genesMatrix, classes, datasetName = "Example",
                    trainParams = TrainParams(fisherDiscriminant, FALSE, TRUE),
                    predictParams = PredictParams(function(){}, FALSE, getClasses = function(result) result),
                    resubstituteParams = ResubstituteParams(nFeatures = seq(10, 100, 10),
                                                            performanceType = "balanced", better = "lower"))
  }
}
