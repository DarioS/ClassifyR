\name{kTSPclassifier}
\alias{kTSPclassifier}
\alias{kTSPclassifier,matrix-method}
\alias{kTSPclassifier,DataFrame-method}
\alias{kTSPclassifier,MultiAssayExperiment-method}
\title{Classification Using k Pairs of Features With Relative Differences Between Classes}
\description{
  Each pair of features votes for a class based on whether the value of one feature is less than
  the other feature.
}
\usage{
  \S4method{kTSPclassifier}{matrix}(measurements, classes, test, featurePairs = NULL, ...)
  \S4method{kTSPclassifier}{DataFrame}(measurements, classes, test, featurePairs = NULL,
                   weighted = c("unweighted", "weighted", "both"),
                   minDifference = 0, returnType = c("both", "class", "score"), verbose = 3)
  \S4method{kTSPclassifier}{MultiAssayExperiment}(measurements, test, target = names(measurements)[1],
                   featurePairs = NULL, ...)
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
  \item{featurePairs}{An object of class as \code{\link{Pairs}} containing the pairs of features to determine
                      whether the inequality of the first feature being less than the second feature holds, 
                      indiciating evidence for the second level of the \code{classes} factor.}              
  \item{target}{If \code{measurements} is a \code{MultiAssayExperiment}, the name of the
                 data table to be used. \code{"clinical"} is also a valid value and specifies that
                 integer variables from the clinical data table will be used.}              
  \item{...}{Unused variables by the methods for a \code{matrix} or a \code{MultiAssayExperiment} passed to the
             \code{DataFrame} method which does the classification.}
  \item{weighted}{Default: \code{"both"}. Either \code{"both"}, \code{"unweighted"} or \code{"weighted"}. 
                  In weighted mode, the difference in densities is summed over all features.
                  If unweighted mode, each feature's vote is worth the same. Both can be calculated simultaneously.}
  \item{minDifference}{Default: 0. The minimum difference in densities for a feature to be allowed
                       to vote. Can be a vector of cutoffs. If no features for a particular sample
                       have a difference large enough, the class predicted is simply the largest class.}
  \item{returnType}{Default: \code{"both"}. Either \code{"class"}, \code{"score"} or \code{"both"}.
                    Sets the return value from the prediction to either a vector of class labels,
                    score for a sample belonging to the second class, as determined by the factor levels,
                    or both labels and scores in a \code{data.frame}.}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  If \code{weighted} is \code{TRUE}, then a sample's predicted class is based on the sum of differences of 
  measurements for each feature pair. Otherwise, when \code{weighted} is \code{FALSE}, each pair of features
  has an equal vote, the predicted class is the one with the most votes. If the voting is tied, the the class
  with the most samples in the training set is voted for.
  
  Because this method compares different features, they need to have comparable measurements. For example,
  RNA-seq counts would be unsuitable since these depend on the length of a feature, whereas
  F.P.K.M. values would be suitable.
  
  The \code{featurePairs} to use is recommended to be determined in conjunction with \code{pairsDifferencesSelection}.
}
\value{
  A vector or list of class prediction information, as long as the number of samples in the test data,
  or lists of such information, if a variety of predictions is generated.
}
\author{Dario Strbenac}
\seealso{\code{\link{pairsDifferencesSelection}} for a function which could be used to do feature selection before
  the k-TSP classifier is run.}
\examples{
  # Difference in differences for features A and C between classes.                                           
  measurements <- matrix(c(9.9, 10.5, 10.1, 10.9, 11.0, 6.6, 7.7, 7.0, 8.1, 6.5,
                           8.5, 10.5, 12.5, 10.5, 9.5, 8.5, 10.5, 12.5, 10.5, 9.5,
                           6.6, 7.7, 7.0, 8.1, 6.5, 11.2, 11.0, 11.1, 11.4, 12.0,
                           8.1, 10.6, 7.4, 7.1, 10.4, 6.1, 7.3, 2.7, 11.0, 9.1,
                           round(rnorm(60, 8, 1), 1)), ncol = 10, byrow = TRUE)
  classes <- factor(rep(c("Good", "Poor"), each = 5))
                         
  rownames(measurements) <- LETTERS[1:10]
  colnames(measurements) <- names(classes) <- paste("Patient", 1:10)
  
  trainIndex <- c(1:4, 6:9)
  trainMatrix <- measurements[, trainIndex]
  testMatrix <- measurements[, c(5, 10)]
  
  featurePairs <- Pairs('A', 'C') # Could be selected by pairsDifferencesSelection function.
  kTSPclassifier(trainMatrix, classes[trainIndex], testMatrix, featurePairs)
}
