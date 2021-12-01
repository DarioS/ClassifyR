\name{CrossValParams}
\docType{class}
\alias{CrossValParams}
\alias{CrossValParams-class}
\alias{CrossValParams,missing,missing,missing,missing,missing,missing,missing,missing-method}
\alias{CrossValParams,characterOrMissing,numericOrMissing,numericOrMissing,numericOrMissing,numericOrMissing,characterOrMissing,BiocParallelParamOrMissing,numericOrMissing-method}

\title{Parameters for Cross-validation Specification}
\description{
  Collects and checks necessary parameters required for cross-validation by \code{\link{runTests}}.
}
\section{Constructor}{
  \describe{
    \item{}{
      \code{CrossValParams()}
      Creates a default \code{CrossValParams} object, specifying that 100 permutations
      and 5 folds of each to be done. A seed for random number generation is set to
      987654321. Features to be selected by resubstitution of the training set for making predictions and
      choosing the best set of features rather than computationally expensive nested cross-validation.
  }
    \item{}{
      \preformatted{CrossValParams(samplesSplits = c("Permute k-Fold", "Permute Percentage Split", "Leave-k-Out", "k-Fold"),
      permutations = 100, percentTest = 25, folds = 5, leave = 2,
      tuneMode = c("Resubstitution", "Nested CV"), parallelParams = bpparam(), seed = 987654321)}
      Creates a \code{ModellingParams} object based on user-specified settings with default values of
      different classification stages geared at homoscedastic numeric data.
	\describe{
            \item{\code{samplesSplits}}{Default: "Permute k-Fold". A character value specifying what kind of sample splitting to do.}
            \item{\code{permutations}}{Default: 100. Number of times to permute the data set before it is split into
            training and test sets. Only relevant if \code{samplesSplits} is either \code{"Permute k-Fold"} or
            \code{"Permute Percentage Split"}.}
            \item{\code{percentTest}}{The percentage of the data set to assign to the test set, with the remainder of the
            samples belonging to the training set. Only relevant if \code{samplesSplits} is  \code{"Permute Percentage Split"}.}
            \item{\code{folds}}{The number of approximately equal-sized folds to partition the samples into. Only relevant if
            \code{samplesSplits} is \code{"Permute k-Fold"} or \code{"k-Fold"}.}
            \item{\code{leave}}{The number of samples to generate all possible combination of and use as the test set.
            Only relevant if \code{samplesSplits} is \code{"Leave-k-Out"}. If set to 1, it is the traditional leave-one-out
            cross-validation, sometimes written as LOOCV.}
            \item{\code{tuneMode}}{Default: Resubstitution. The scheme to use for selecting any tuning parameters.}
            \item{\code{parallelParams}}{An instance of \code{\link{BiocParallelParam}} specifying the kind of parallelisation
            to use. Default is to use one core on Windows operating system computers or two cores
            less than the total number of cores the computer has, if it has four or more cores and
            is an operating system other than Windows, as is the default of \code{\link{bpparam}}.}
            \item{\code{seed}}{Default: 987654321. Random number generator seed. Affects how samples
            are permuted for schemes which use a permutation of samples.}

	         }
    }
  }
}

\author{Dario Strbenac}
\examples{
  CrossValParams() # Default is 100 permutations and 5 folds of each.
  CrossValParams("Leave-k-Out", leave = 2) # Leave-2-out cross-validation.
}