\name{distribution}
\alias{distribution}
\alias{distribution,ClassifyResult-method}
\title{Get Frequencies of Feature Selection and Sample Errors}
\description{There are two modes. For aggregating feature selection results, the function
counts the number of times each feature was selected in all cross validations. For aggregating
classification results, the error rate for each sample is calculated. This is useful in
identifying outlier samples that are difficult to classify.
}
\usage{
  \S4method{distribution}{ClassifyResult}(result, type = c("features", "samples"),
                         summary = c("density", "frequency"), plot = TRUE, xMax = NULL, ...)
}
\arguments{
  \item{result}{An object of class \code{\link{ClassifyResult}}.}
  \item{type}{Whether to calculate sample-wise error rate or the number of times a
              feature was selected.}
  \item{summary}{Whether to plot frequencies or densities. If feature distribution is analysed,
                it will also cause the retured vector to be a decimal representing the percentage.}              
  \item{plot}{Whether to draw a histogram of the aggregation.}
  \item{xMax}{Maximum bin value for histogram to plot.}            
  \item{...}{Further parameters, such as \code{colour} and \code{fill}, passed to
             \code{\link[ggplot2]{geom_histogram}}.}
}
\value{
  If \code{type} is "features", a vector as long as the number of features that were
  chosen at least once containing the number of times the feature was chosen in cross
  validations. If \code{type} is "samples", a vector as long as the number of samples,
  containing the cross validation error rate of the sample.
}
\author{Dario Strbenac}

\examples{
  if(require(curatedOvarianData) && require(sparsediscrim))
  {
    data(TCGA_eset)
    badOutcome <- which(pData(TCGA_eset)[, "vital_status"] == "deceased" & pData(TCGA_eset)[, "days_to_death"] <= 365)
    goodOutcome <- which(pData(TCGA_eset)[, "vital_status"] == "living" & pData(TCGA_eset)[, "days_to_death"] >= 365 * 5)
    TCGA_eset <- TCGA_eset[, c(badOutcome, goodOutcome)]
    classes <- factor(rep(c("Poor", "Good"), c(length(badOutcome), length(goodOutcome))))
    pData(TCGA_eset)[, "class"] <- classes
    result <- runTests(TCGA_eset, resamples = 2, fold = 2)
    sampleDistribution <- distribution(result, "samples", binwidth = 0.1)
    featureDistribution <- distribution(result, "features", binwidth = 1)
    print(head(sampleDistribution))
    print(head(featureDistribution))
  }
}