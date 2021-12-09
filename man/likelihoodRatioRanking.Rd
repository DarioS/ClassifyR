\name{likelihoodRatioRanking}
\alias{likelihoodRatioRanking}
\alias{likelihoodRatioRanking,matrix-method}
\alias{likelihoodRatioRanking,DataFrame-method}
\alias{likelihoodRatioRanking,MultiAssayExperiment-method}
\title{Ranking of Differential Distributions with Likelihood Ratio Statistic}
\description{Ranks features from largest difference of log likelihoods (null hypothesis - alternate hypothesis) to smallest.
}
\usage{
  \S4method{likelihoodRatioRanking}{matrix}(measurements, classes, ...)
  \S4method{likelihoodRatioRanking}{DataFrame}(measurements, classes,
                                alternative = c(location = "different", scale = "different"),
                                ..., verbose = 3)
  \S4method{likelihoodRatioRanking}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
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
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 numeric variables from the clinical data table will be used.}
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method or options which are accepted by
             the function \code{\link{getLocationsAndScales}}.}
  \item{alternative}{Default: \code{c("different", "different")}. A vector of length 2.
                    The first element specifies the location of the alternate hypothesis.
                    The second element specifies the scale of the alternate hypothesis.
                    Valid values in each element are "same" or "different".}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  Likelihood ratio test of null hypothesis that the location and scale are the same for
  both groups, and an alternate hypothesis that is specified by parameters. The location and scale
  of features is calculated by \code{\link{getLocationsAndScales}}. The distribution fitted to the data
  is the normal distribution.
}
\value{
  A vector or data frame (if \code{MultiAssayExperiment} input) of features, from the most promising features in the first position to the least promising feature in the last position.
}
\author{Dario Strbenac}

\examples{
  # First 20 features have bimodal distribution for Poor class.
  # Other 80 features have normal distribution for both classes.

  genesMatrix <- sapply(1:25, function(sample)
                              {
                                randomMeans <- sample(c(8, 12), 20, replace = TRUE)
                                c(rnorm(20, randomMeans, 1), rnorm(80, 10, 1))
                              }
                       )
  genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample) rnorm(100, 10, 1)))
  rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
  classes <- factor(rep(c("Poor", "Good"), each = 25))

  ranked <- likelihoodRatioRanking(genesMatrix, classes)
  head(ranked)
}
