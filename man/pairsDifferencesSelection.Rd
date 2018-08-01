\name{pairsDifferencesSelection}
\alias{pairsDifferencesSelection}
\alias{pairsDifferencesSelection,matrix-method}
\alias{pairsDifferencesSelection,DataFrame-method}
\alias{pairsDifferencesSelection,MultiAssayExperiment-method}
\title{Selection of Pairs of Features that are Different Between Classes}
\description{Ranks pre-specified pairs of features by the largest difference of the sum of
             measurement differences over all samples within a class and chooses the
             pairs of features which have the best resubstitution performance.
}
\usage{
  \S4method{pairsDifferencesSelection}{matrix}(measurements, classes, featurePairs = NULL, ...)
  \S4method{pairsDifferencesSelection}{DataFrame}(measurements, classes, featurePairs = NULL,
                  datasetName, trainParams, predictParams, resubstituteParams,
                  selectionName = "Pairs Differences", verbose = 3)
  \S4method{pairsDifferencesSelection}{MultiAssayExperiment}(measurements, target = names(measurements)[1], featurePairs = NULL, ...)
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
  \item{featurePairs}{An S4 object of type Pairs containing feature identifiers to calculate
                      the sum of differences within each class for.}
  \item{target}{If \code{measurements} is a \code{MultiAssayExperiment}, the name of the
                 data table to be used.}
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method.}
  \item{datasetName}{A name for the data set used. Stored in the result.}
  \item{trainParams}{A container of class \code{\link{TrainParams}} describing the
    classifier to use for training.}
  \item{predictParams}{A container of class \code{\link{PredictParams}} describing how
    prediction is to be done.}
  \item{resubstituteParams}{An object of class \code{\link{ResubstituteParams}}
    describing the performance measure to consider and the numbers of
    top features to try for resubstitution classification.}                               
  \item{selectionName}{A name to identify this selection method by. Stored in the result.}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
    This function only prints progress messages if the value is 3.}
}
\details{
  Instead of considering whether one feature in a pair of features is consistently lower or higher
  than the other in the pair, this method takes the sum of differences across all samples
  within a class, to prevent ties in the ranking of pairs of features.
}
\value{
  An object of class \code{\link{SelectResult}} or a list of such objects, if the classifier which was
  used for determining the specified performance metric made a number of prediction varieties.
}
\author{Dario Strbenac}
\seealso{
\code{\link{kTSPclassifier}} for a classifier which makes use of the pairs of selected features in classification.
}
\references{
  Simple decision rules for classifying human cancers from gene expression profiles,
  Aik C Tan, Daniel Q Naiman, Lei Xu, Raimond L. Winslow and Donald Geman, 2005,
  \emph{Bioinformatics}, Volume 21 Issue 20, \url{https://academic.oup.com/bioinformatics/article/21/20/3896/203010}.
}
\examples{
  featurePairs <- Pairs(c('A', 'A'), c('B', 'C'))
                           
  # Difference in differences for features A and C between classes.                                           
  measurements <- matrix(c(9.9, 10.5, 10.1, 10.9, 11.0, 6.6, 7.7, 7.0, 8.1, 6.5,
                           8.5, 10.5, 12.5, 10.5, 9.5, 8.5, 10.5, 12.5, 10.5, 9.5,
                           6.6, 7.7, 7.0, 8.1, 6.5, 11.2, 11.0, 11.1, 11.4, 12.0,
                           8.1, 10.6, 7.4, 7.1, 10.4, 6.1, 7.3, 2.7, 11.0, 9.1,
                           round(rnorm(60, 8, 1), 1)), ncol = 10, byrow = TRUE)
  classes <- factor(rep(c("Good", "Poor"), each = 5))
                         
  rownames(measurements) <- LETTERS[1:10]
  colnames(measurements) <- names(classes) <- paste("Patient", 1:10)

  # The features are pairs and there are only two in this example.
  resubstituteParams <- ResubstituteParams(nFeatures = 1:2,
                                performanceType = "balanced error", better = "lower")
  
  predictParams <- PredictParams(NULL, getClasses = function(result) result)
  pairsDifferencesSelection(measurements, classes, featurePairs = featurePairs,
                            datasetName = "Example",
                            trainParams = TrainParams(kTSPclassifier),
                            predictParams = predictParams,
                            resubstituteParams = resubstituteParams)
}