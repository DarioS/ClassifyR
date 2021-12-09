\name{KolmogorovSmirnovRanking}
\alias{KolmogorovSmirnovRanking}
\alias{KolmogorovSmirnovRanking,matrix-method}
\alias{KolmogorovSmirnovRanking,DataFrame-method}
\alias{KolmogorovSmirnovRanking,MultiAssayExperiment-method}
\title{Ranking of Differential Distributions with Kolmogorov-Smirnov Distance}
\description{Ranks features from largest Kolmogorov-Smirnov distance to smallest.
}
\usage{
  \S4method{KolmogorovSmirnovRanking}{matrix}(measurements, classes, ...)
  \S4method{KolmogorovSmirnovRanking}{DataFrame}(measurements, classes, ...,
                                   verbose = 3)
  \S4method{KolmogorovSmirnovRanking}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
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
             are passed into and used by the \code{DataFrame} method or options which are accepted by the function
             \code{\link{ks.test}}.}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
               This function only prints progress messages if the value is 3.}
}
\value{
  A vector or data frame (if \code{MultiAssayExperiment} input) of features, from the most promising features in the first position to the least promising feature in the last position.
}
\author{Dario Strbenac}

\examples{
  # First 20 features have bimodal distribution for Poor class.
  # Other 80 features have normal distribution for both classes.
  set.seed(1984)
  genesMatrix <- sapply(1:25, function(sample)
                              {
                                randomMeans <- sample(c(8, 12), 20, replace = TRUE)
                                c(rnorm(20, randomMeans, 1), rnorm(80, 10, 1))
                              }
                       )
  genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample) rnorm(100, 10, 1)))
  rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
  classes <- factor(rep(c("Poor", "Good"), each = 25))

  ranked <- KolmogorovSmirnovRanking(genesMatrix, classes)
  head(ranked)                           
}
