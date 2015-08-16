\name{previousSelection}
\alias{previousSelection}
\alias{previousSelection,matrix-method}
\alias{previousSelection,ExpressionSet-method}
\title{Automated Selection of Previously Selected Features}
\description{Uses the feature selection of the same cross-validation iteration of a previous classification
for the current classification task.
}
\usage{
  \S4method{previousSelection}{matrix}(expression, classes, ...)
  \S4method{previousSelection}{ExpressionSet}(expression, datasetName, classifyResult,
                                              minimumOverlapPercent = 80,
                                              selectionName = "Previous Selection", .iteration, verbose = 3)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
    the training data. For a matrix, the rows are features, and the columns are samples.}
  \item{classes}{A vector of class labels.}
  \item{...}{For the \code{\link{matrix}} method, variables passed to the \code{\link{ExpressionSet}} method.}  
  \item{datasetName}{A name for the dataset used. Stored in the result.}
  \item{classifyResult}{An existing classification result from which to take the feature selections from.}
  \item{minimumOverlapPercent}{If at least this many selected features can't be identified in the
                               current dataset, then the selection stops with an error.}
  \item{selectionName}{A name to identify this selection method by. Stored in the result.}
  \item{.iteration}{Not to be set by the user.}
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
    This function only prints progress messages if the value is 3.}
}
\value{
  An object of class \code{\link{SelectResult}}.
}
\author{Dario Strbenac}

\examples{
  if(require(sparsediscrim))
  {
    # Genes 76 to 100 have differential expression.
    genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
                                 c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
    rownames(genesMatrix) <- paste("Gene", 1:100)                                 
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    resubstitute <- ResubstituteParams(nFeatures = seq(10, 100, 10), 
                                       performanceType = "err", better = "lower")
    result <- runTests(genesMatrix, classes, "Ovarian Cancer", "Differential Expression",
                       resamples = 2, fold = 2,
                       params = list(SelectParams(limmaSelection, resubstituteParams = resubstitute),
                       TrainParams(dlda, TRUE, FALSE),
                       PredictParams(predict, TRUE, getClasses = function(result) result[["class"]])))
                       
    # Genes 74 to 98 have differential expression in new dataset.
    newDataset <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
    newDataset <- cbind(newDataset, rbind(sapply(1:25, function(sample) rnorm(73, 9, 2)),
                                           sapply(1:25, function(sample) rnorm(25, 14, 2)),
                                           sapply(1:25, function(sample) rnorm(2, 14, 2))))
                                           
    newerResult <- runTests(newDataset, classes, "Ovarian Cancer Updated", "Differential Expression",
                       resamples = 2, fold = 2,
                       params = list(SelectParams(previousSelection, intermediate = ".iteration",
                                    classifyResult = result),
                       TrainParams(dlda, TRUE, FALSE),
                       PredictParams(predict, TRUE, getClasses = function(result) result[["class"]])))                                           
  }  
}
