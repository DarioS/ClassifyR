\name{fisherDiscriminant}
\alias{fisherDiscriminant}
\alias{fisherDiscriminant,matrix-method}
\alias{fisherDiscriminant,ExpressionSet-method}
\title{Classification Using Fisher's LDA}
\description{Finds the decision boundary using the training set, and gives predictions
for the test set.
}
\usage{
  \S4method{fisherDiscriminant}{matrix}(expression, classes, test, verbose = 3)
  \S4method{fisherDiscriminant}{ExpressionSet}(expression, test, verbose = 3)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    the training data. For a matrix, the rows are features, and the columns
                    are samples.}
  \item{classes}{A vector of class labels.}
  \item{test}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    the test data.}
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  Unlike ordinary LDA, Fisher's version does not have assumptions about the normality
  of the features.
}
\value{
  A vector of class predictions, as long as the number of samples in the test data.
}
\author{Dario Strbenac}

\examples{
  trainMatrix <- matrix(rnorm(1000, 8, 2), ncol = 10)
  trainMatrix[1:30, 1:5] <- trainMatrix[1:30, 1:5] + 5 # Make first 30 genes D.E.
  testMatrix <- matrix(rnorm(1000, 8, 2), ncol = 10)
  testMatrix[1:30, 6:10] <- testMatrix[1:30, 6:10] + 5 # Make first 30 genes D.E.
  classes <- factor(rep(c("Poor", "Good"), each = 5))
  fisherDiscriminant(trainMatrix, classes, testMatrix)
}
