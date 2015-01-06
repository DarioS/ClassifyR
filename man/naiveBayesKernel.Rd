\name{naiveBayesKernel}
\alias{naiveBayesKernel}
\alias{naiveBayesKernel,matrix-method}
\alias{naiveBayesKernel,ExpressionSet-method}
\title{Classification Using A Bayes Classifier with Kernel Density Estimates}
\description{
  Kernel density estimates are fitted to the training data and a naive Bayes classifier
  is used to classify samples in the test data. 
}
\usage{
  \S4method{naiveBayesKernel}{matrix}(expression, classes, ...)
  \S4method{naiveBayesKernel}{ExpressionSet}(expression, test,
            weighted = c("both", "unweighted", "weighted"), minDifference = 0,
            returnType = c("labels", "scores", "both"), verbose = 3)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    the training data. For a matrix, the rows are features, and the columns
                    are samples.}
  \item{classes}{A vector of class labels.}
  \item{...}{Unused variables from the \code{\link{matrix}} method passed to the
             \code{\link{ExpressionSet}} method.}            
  \item{test}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    the test data.}
  \item{weighted}{In weighted mode, the difference in densities is summed over all features.
                  If unweighted mode, each features's vote is worth the same. To save
                  computational time, both can be calculated simultaneously.}
  \item{minDifference}{The minimum difference in densities for a feature to be allowed
                       to vote. Can be a vector of cutoffs.}
  \item{returnType}{Either \code{"labels"}, \code{"scores"}, or \code{"both"}. Sets the return value
                    from the prediction to either a vector of class labels, score for a sample belonging
                    to the second class, as determined by the factor levels, or both labels and scores
                    in a \code{\link{data.frame}}.}     
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  If \code{weighted} is \code{TRUE}, then a sample's predicted class is the class with
  the largest sum of density estimates, scaled for the number of samples in
  the training data of each class. Otherwise, when \code{weighted} is \code{FALSE},
  each feature has an equal vote, and votes for the class with the largest density,
  scaled for class sizes in the training set.
  
  For unweighted classification, the score is the number of features which voted for the second class,
  divided by the features which voted. For weighted classification, the score is the sum of scaled
  differences in densities.
}
\value{
  A vector or list of class prediction information, as long as the number of samples in the test data,
  or lists of such information, if both weighted and unweighted voting or a range of \code{minDifference}
  values was provided.
}
\author{Dario Strbenac, John Ormerod}
\examples{
  trainMatrix <- matrix(rnorm(1000, 8, 2), ncol = 10)
  trainMatrix[1:30, 1:5] <- trainMatrix[1:30, 1:5] + 5 # Make first 30 genes D.E.
  testMatrix <- matrix(rnorm(1000, 8, 2), ncol = 10)
  testMatrix[1:30, 6:10] <- testMatrix[1:30, 6:10] + 5 # Make first 30 genes D.E.
  classes <- factor(rep(c("Poor", "Good"), each = 5))
  naiveBayesKernel(trainMatrix, classes, testMatrix)
}
