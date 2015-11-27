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
  \S4method{distribution}{ClassifyResult}(result, dataType = c("features", "samples"),
                   plotType = c("density", "histogram"), summaryType = c("percentage", "count"),
                   plot = TRUE, xMax = NULL, xLabel = "Percentage of Cross-validations",
                   yLabel = "Density", title = "Distribution of Feature Selections",
                   fontSizes = c(24, 16, 12), ...)
}
\arguments{
  \item{result}{An object of class \code{\link{ClassifyResult}}.}
  \item{dataType}{Whether to calculate sample-wise error rate or the number of times a
              feature was selected.}
  \item{plotType}{Whether to draw a probability density curve or a histogram.}
  \item{summaryType}{Whether to summarise the feature selections as a percentage or count.}
  \item{plot}{Whether to draw a plot of the frequency of selection or error rate.}
  \item{xMax}{Maximum data value to show in plot.}
  \item{xLabel}{The label for the x-axis of the plot.}
  \item{yLabel}{The label for the y-axis of the plot.}
  \item{title}{An overall title for the plot.}
  \item{fontSizes}{A vector of length 3. The first number is the size of the title.
                   The second number is the size of the axes titles. The third number is
                   the size of the axes values.}
  \item{...}{Further parameters, such as \code{colour} and \code{fill}, passed to
             \code{\link[ggplot2]{geom_histogram}} or \code{\link[ggplot2]{stat_density}},
             depending on the value of \code{plotType}.}
}
\value{
  If \code{type} is "features", a vector as long as the number of features that were
  chosen at least once containing the number of times the feature was chosen in cross
  validations or the percentage of times chosen. If \code{type} is "samples",
  a vector as long as the number of samples, containing the cross-validation error rate
  of the sample. If \code{plot} is TRUE, then a plot is also made on the current
  graphics device.
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
    result <- runTests(TCGA_eset, "Ovarian Cancer", "Differential Expression", resamples = 2, fold = 2)
    sampleDistribution <- distribution(result, "samples", xLabel = "Sample Error Rate",
                                       title = "Distribution of Error Rates")
    featureDistribution <- distribution(result, "features", summaryType = "count", plotType = "histogram",
                                        xLabel = "Number of Cross-validations", yLabel = "Count",
                                        binwidth = 1)
    print(head(sampleDistribution))
    print(head(featureDistribution))
  }
}