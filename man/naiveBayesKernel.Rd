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
  \S4method{naiveBayesKernel}{ExpressionSet}(expression, test, densityFunction = density,
            densityParameters = list(bw = "SJ", n = 1024, from = expression(min(featureValues)),
                                                          to = expression(max(featureValues))),
            weighted = c("both", "unweighted", "weighted"),
            weight = c("both", "height difference", "crossover distance"),
            minDifference = 0, tolerance = 0.01, returnType = c("label", "score", "both"), verbose = 3)
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
  \item{densityFunction}{A function which will return a probability density, which is essentially
                         a list with x and y coordinates.}
  \item{densityParameters}{A list of options for \code{densityFunction}.}                       
  \item{weighted}{In weighted mode, the difference in densities is summed over all features.
                  If unweighted mode, each feature's vote is worth the same. To save
                  computational time, both can be calculated simultaneously.}
  \item{weight}{The type of weight to calculate. For \code{"height difference"}, the weight of each prediction
                is equal to the verical distance between two densities, for a particular value of x. For
                \code{"height difference"}, the x positions where two densities cross is firstly calculated.
                The predicted class is the class with the highest density at the particular value of x and
                the weight is the distance of x from the nearest density crossover point.}               
  \item{minDifference}{The minimum difference in densities for a feature to be allowed
                       to vote. Can be a vector of cutoffs. If no features for a particular sample
                       have a difference large enough, the class predicted is simply the largest class.}
  \item{tolerance}{Only relevant when \code{weight} is \code{"crossover distance"}. Absolute differences in
                   the y values of the two densities of this magnitude or smaller cause the densities at the
                   corresponding x values to be considered as overlapping.}
  \item{returnType}{Either \code{"label"}, \code{"score"}, or \code{"both"}. Sets the return value
                    from the prediction to either a vector of class labels, score for a sample belonging
                    to the second class, as determined by the factor levels, or both labels and scores
                    in a \code{\link{data.frame}}.}     
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  If \code{weighted} is \code{TRUE}, then a sample's predicted class is the class with
  the largest sum of weights, scaled for the number of samples in
  the training data of each class. Otherwise, when \code{weighted} is \code{FALSE},
  each feature has an equal vote, and votes for the class with the largest weight,
  scaled for class sizes in the training set.
  
  For unweighted classification, the score is the number of features which voted for the second class,
  divided by the features which voted. For weighted classification, the score is the sum of scaled
  differences of second class densities to first class densities.
  
  The variable name of each feature's measurements in the iteration over all features is \code{featureValues}. 
  This is important to know if each feature's measurements need to be referred to in the specification of
  \code{densityParameters}, such as for specifying the range of x values of the density function to be computed.
  
  If \code{weight} is \code{"crossover distance"}, the crossover points are computed by considering the 
  distance between y values of the two densities at every x value. If the y values are sufficiently close, the
  corresponding x values added to a candidate list. Consecutive x values are grouped, and the x value in each group
  that has with the smallest distance is chosen as the representative location of the crossover point. Only y values   that are \code{tolerance} or greater are considered in this first stage. If no crossover points are found,
  the y values below \code{tolerance} are considered, except for those at the leftmost or rightmost region of
  the range of the densities. This is necessary when the densities are completely separated.
}
\value{
  A vector or list of class prediction information, as long as the number of samples in the test data,
  or lists of such information, if a variety of predictions is generated.
}
\author{Dario Strbenac, John Ormerod}
\examples{
  trainMatrix <- matrix(rnorm(1000, 8, 2), ncol = 10)
  trainMatrix[1:30, 1:5] <- trainMatrix[1:30, 1:5] + 5 # Make first 30 genes D.E.
  testMatrix <- matrix(rnorm(1000, 8, 2), ncol = 10)
  testMatrix[1:30, 6:10] <- testMatrix[1:30, 6:10] + 5 # Make first 30 genes D.E.
  classes <- factor(rep(c("Poor", "Good"), each = 5))
  # Expected: Good Good Good Good Good Poor Poor Poor Poor Poor
  naiveBayesKernel(trainMatrix, classes, testMatrix)
}
