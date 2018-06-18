\name{featureSetSummary}
\alias{featureSetSummary}
\alias{featureSetSummary,matrix-method}
\alias{featureSetSummary,DataFrame-method}
\alias{featureSetSummary,MultiAssayExperiment-method}
\title{Transform a Table of Feature Abundances into a Table of Feature Set Abundances.}
\description{Represents a feature set by the mean or median feature measurement of a feature set for
             all features belonging to a feature set.
             }
\usage{
  \S4method{featureSetSummary}{matrix}(measurements, location = c("median", "mean"),
                   featureSets, minimumOverlapPercent = 80, verbose = 3)
  \S4method{featureSetSummary}{DataFrame}(measurements, location = c("median", "mean"),
                                          featureSets, minimumOverlapPercent = 80, verbose = 3)
  \S4method{featureSetSummary}{MultiAssayExperiment}(measurements, target = NULL, location = c("median", "mean"),
                   featureSets, minimumOverlapPercent = 80, verbose = 3)
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}} or \code{\link{DataFrame}} containing
                      the training data. For a \code{matrix}, the rows are features, and
                      the columns are samples.}
  \item{target}{If the input is a \code{\link{MultiAssayExperiment}}, this specifies which
                data set will be transformed. Can either be an integer or a character string
                specifying the name of the table. Must have length 1.}
  \item{location}{Default: The median. The type of location to summarise a set of features
                  belonging to a feature set by.}
  \item{featureSets}{An object of type \code{\link{FeatureSetCollection}} which defines the
                     feature sets.}
  \item{minimumOverlapPercent}{The minimum percentage of overlapping features between the data set
                               and a feature set defined in \code{featureSets} for that
                               feature set to not be discarded from the anaylsis.}                     
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  This feature transformation method is unusual because the mean or median feature of a feature set for
  one sample may be different to another sample, whereas most other feature transformation methods do
  not result in different features being compared between samples during classification.
}
\value{
  The same class of variable as the input variable \code{measurements} is, with the individual features
  summarised to feature sets. The number of samples remains unchanged, so only one dimension of
  \code{measurements} is altered.
}

\author{Dario Strbenac}
\references{
  Network-based biomarkers enhance classical approaches to prognostic gene expression signatures,
  Rebecca L Barter, Sarah-Jane Schramm, Graham J Mann and Yee Hwa Yang, 2014,
  \emph{BMC Systems Biology}, Volume 8 Supplement 4 Article S5, \url{https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-8-S4-S5}.
}
\examples{
  sets <- list(Adhesion = c("Gene 1", "Gene 2", "Gene 3"),
               `Cell Cycle` = c("Gene 8", "Gene 9", "Gene 10"))
  featureSets <- FeatureSetCollection(sets)
  
  # Adhesion genes have a median gene difference between classes.
  genesMatrix <- matrix(c(rnorm(5, 9, 0.3), rnorm(5, 7, 0.3), rnorm(5, 8, 0.3),
                        rnorm(5, 6, 0.3), rnorm(10, 7, 0.3), rnorm(70, 5, 0.1)),
                        ncol = 10, byrow = TRUE)
  rownames(genesMatrix) <- paste("Gene", 1:10)
  colnames(genesMatrix) <- paste("Patient", 1:10)
  classes <- factor(rep(c("Poor", "Good"), each = 5)) # But not used for transformation.
  
  featureSetSummary(genesMatrix, featureSets = featureSets)
}
