\name{calcPerformance}
\alias{calcPerformance}
\alias{calcExternalPerformance}
\alias{calcCVperformance}
\alias{calcExternalPerformance,factor,factor-method}
\alias{calcCVperformance,ClassifyResult-method}
\title{Add Performance Calculations to a ClassifyResult Object or Calculate for a Pair of Factor Vectors}
\description{If \code{calcExternalPerformance} is used, such as when having a vector of known classes and a vector of predicted classes determined outside of the ClassifyR package, a single metric value is calculated. If \code{calcCVperformance} is used, annotates the results of calling \code{\link{runTests}} with one of the user-specified performance measures.
}
\usage{
  \S4method{calcExternalPerformance}{factor,factor}(actualClasses, predictedClasses,
                              performanceType = c("error", "balanced error", "average accuracy",
                                       "micro precision", "micro recall",
                                       "micro F1", "macro precision",
                                       "macro recall", "macro F1"))
  \S4method{calcCVperformance}{ClassifyResult}(result, performanceType = c("error", "accuracy", "balanced error", "balanced accuracy",
                                               "sample error", "sample accuracy",
                                              "micro precision", "micro recall",
                                             "micro F1", "macro precision",
                                             "macro recall", "macro F1"))
}
\arguments{
  \item{result}{An object of class \code{\link{ClassifyResult}}.}
  \item{performanceType}{A character vector of length 1. Default: \code{"balanced error"}.\cr
                         Must be one of the following options:\cr
                         \itemize{
                         \item{\code{"error"}: Ordinary error rate.}
                         \item{\code{"accuracy"}: Ordinary accuracy.}
                         \item{\code{"balanced error"}: Balanced error rate.}
                         \item{\code{"balanced accuracy"}: Balanced accuracy.}
                         \item{\code{"sample error"}: Error rate for each sample in the data set.}
                         \item{\code{"sample accuracy"}: Accuracy for each sample in the data set.}

                         \item{\code{"micro precision"}: Sum of the number of correct predictions in each class,
                                                   divided by the sum of number of samples
                                                   in each class.}
                         \item{\code{"micro recall"}: Sum of the number of correct predictions in each class,
                                                divided by the sum of number of samples predicted
                                                as belonging to each class.}
                         \item{\code{"micro F1"}: F1 score obtained by calculating the harmonic mean of
                                            micro precision and micro recall.}
                         \item{\code{"macro precision"}: Sum of the ratios of the number of correct predictions
                                                   in each class to the number of samples in each class,
                                                   divided by the number of classes.}
                         \item{\code{"macro recall"}: Sum of the ratios of the number of correct predictions
                                                in each class to the number of samples predicted to be
                                                in each class, divided by the number of classes.}
                         \item{\code{"macro F1"}: F1 score obtained by calculating the harmonic mean of
                                            macro precision and macro recall.}
                        }
                        }
  \item{actualClasses}{A factor vector specifying each sample's correct class.}
  \item{predictedClasses}{A factor vector of the same length as \code{actualClasses}
                          specifying each sample's predicted class.}
}
\details{
  All metrics are suitable for evaluating classification scenarios with more than two classes
  and are reimplementations of those available from \href{https://software.intel.com/en-us/daal-programming-guide-details-40}{Intel DAAL}.

  If \code{\link{runTests}} was run in resampling mode, one performance measure is produced
  for every resampling. If the leave-k-out mode was used, then the predictions are
  concatenated, and one performance measure is calculated for all classifications.
  
  \code{"balanced error"} calculates the balanced error rate and is better suited to class-imbalanced data sets
  than the ordinary error rate specified by {\code{"error"}}. \code{"sample error"} calculates the error rate
  of each sample individually. This may help to identify which samples are contributing the most to the
  overall error rate and check them for confounding factors. Precision, recall and F1 score have micro and macro
  summary versions. The macro versions are preferable because the metric will not have a good score if there is
  substantial class imbalance and the classifier predicts all samples as belonging to the majority class.
}
\value{
  An updated \code{\linkS4class{ClassifyResult}} object, with new information in the \code{performance} slot.
}
\examples{
  predictTable <- data.frame(sample = 1:10,
                             class = factor(sample(LETTERS[1:2], 50, replace = TRUE)))
  actual <- factor(sample(LETTERS[1:2], 10, replace = TRUE))                             
  result <- ClassifyResult("Example", "Differential Expression", "A Selection",
                           paste("A", 1:10, sep = ''), paste("Gene", 1:50, sep = ''),
                           50, list(1:50, 1:50), list(1:5, 6:15),
                           list(predictTable), actual, list("leave", 2))
  result <- calcCVperformance(result, "balanced error") 
  performance(result)
}
\author{Dario Strbenac}
