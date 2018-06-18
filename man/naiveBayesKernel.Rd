\name{naiveBayesKernel}
\alias{naiveBayesKernel}
\alias{naiveBayesKernel,matrix-method}
\alias{naiveBayesKernel,DataFrame-method}
\alias{naiveBayesKernel,MultiAssayExperiment-method}
\title{Classification Using A Bayes Classifier with Kernel Density Estimates}
\description{
  Kernel density estimates are fitted to the training data and a naive Bayes classifier
  is used to classify samples in the test data.
}
\usage{
  \S4method{naiveBayesKernel}{matrix}(measurements, classes, test, ...)
  \S4method{naiveBayesKernel}{DataFrame}(measurements, classes, test,
                 densityFunction = density,
                 densityParameters = list(bw = "nrd0", n = 1024,
                                                 from = expression(min(featureValues)),
                                              to = expression(max(featureValues))),
                   weighted = c("both", "unweighted", "weighted"),
                   weight = c("both", "height difference", "crossover distance"),
                   minDifference = 0, returnType = c("class", "score", "both"), verbose = 3)
  \S4method{naiveBayesKernel}{MultiAssayExperiment}(measurements, test, targets = names(measurements), ...)  
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
                       \code{\link{MultiAssayExperiment}} containing the training data.
                       For a \code{matrix}, the rows are features, and the columns are samples.}
  \item{classes}{Either a vector of class labels of class \code{\link{factor}} of the same length
                 as the number of samples in \code{measurements} or if the measurements are
                 of class \code{DataFrame} a character vector of length 1 containing the
                 column name in \code{measurement} is also permitted. Not used if \code{measurements}
                 is a \code{MultiAssayExperiment} object.}
  \item{test}{An object of the same class as \code{measurements} with no samples in common with
              \code{measurements} and the same number of features as it.}
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 integer variables from the clinical data table will be used.}              
  \item{...}{Unused variables by the three top-level methods passed to the internal method
             which does the classification.}
  \item{densityFunction}{Default: \code{\link{density}}. A function which will return a probability density,
                         which is essentially a list with x and y coordinates.}
  \item{densityParameters}{A list of options for \code{densityFunction}.
                          Default: \code{list(bw = "nrd0", n = 1024,
                                              from = expression(min(featureValues)),
                                              to = expression(max(featureValues))}.}                       
  \item{weighted}{Default: \code{"both"}. Either \code{"both"}, \code{"unweighted"} or \code{"weighted"}. 
                  In weighted mode, the difference in densities is summed over all features.
                  If unweighted mode, each feature's vote is worth the same. Both can be calculated simultaneously.}
  \item{weight}{Default: \code{"both"}. Either \code{"both"}, \code{"height difference"},
                or \code{"crossover distance"}. The type of weight to calculate.
                For \code{"height difference"}, the weight of each prediction is equal to the
                vertical distance between two densities, for a particular value of x. For
                \code{"crossover distance"}, the x positions where two densities cross is firstly calculated.
                The predicted class is the class with the highest density at the particular value of x and
                the weight is the distance of x from the nearest density crossover point.              
  \item{minDifference}{Default: 0. The minimum difference in densities for a feature to be allowed
                       to vote. Can be a vector of cutoffs. If no features for a particular sample
                       have a difference large enough, the class predicted is simply the largest class.}
  \item{returnType}{Default: \code{"class"}. Either \code{"class"}, \code{"score"} or \code{"both"}.
                    Sets the return value from the prediction to either a vector of class labels,
                    score for a sample belonging to the second class, as determined by the factor levels,
                    or both labels and scores in a \code{data.frame}.}     
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  If \code{weighted} is \code{TRUE}, then a sample's predicted class is the class with
  the largest sum of weights, each scaled for the number of samples in
  the training data of each class. Otherwise, when \code{weighted} is \code{FALSE},
  each feature has an equal vote, and votes for the class with the largest weight,
  scaled for class sizes in the training set.
  
  The variable name of each feature's measurements in the iteration over all features is \code{featureValues}. 
  This is important to know if each feature's measurements need to be referred to in the specification of
  \code{densityParameters}, such as for specifying the range of x values of the density function to be computed.
  For example, see the default value of \code{densityParameters} above.
  
  If \code{weight} is \code{"crossover distance"}, the crossover points are computed by considering the 
  distance between y values of the two densities at every x value. x values for which the sign of the difference
  changes compared to the difference of the closest lower value of x are used as the crossover points.
}
\value{
  A vector or list of class prediction information, as long as the number of samples in the test data,
  or lists of such information, if a variety of predictions is generated.
}
\author{Dario Strbenac, John Ormerod}
\examples{
  trainMatrix <- matrix(rnorm(1000, 8, 2), ncol = 10)
  classes <- factor(rep(c("Poor", "Good"), each = 5))
  
  # Make first 30 genes increased in value for poor samples.
  trainMatrix[1:30, 1:5] <- trainMatrix[1:30, 1:5] + 5
  
  testMatrix <- matrix(rnorm(1000, 8, 2), ncol = 10)
  
  # Make first 30 genes increased in value for sixth to tenth samples.
  testMatrix[1:30, 6:10] <- testMatrix[1:30, 6:10] + 5
  
  naiveBayesKernel(trainMatrix, classes, testMatrix)
}
