\name{KolmogorovSmirnovSelection}
\alias{KolmogorovSmirnovSelection}
\alias{KolmogorovSmirnovSelection,matrix-method}
\alias{KolmogorovSmirnovSelection,DataFrame-method}
\alias{KolmogorovSmirnovSelection,MultiAssayExperiment-method}
\title{Selection of Differential Distributions with Kolmogorov-Smirnov Distance}
\description{Ranks features by largest Kolmogorov-Smirnov distance and chooses
             the features which have best resubstitution performance.
}
\usage{
  \S4method{KolmogorovSmirnovSelection}{matrix}(measurements, classes, ...)
  \S4method{KolmogorovSmirnovSelection}{DataFrame}(measurements, classes, datasetName,
                                   trainParams, predictParams, resubstituteParams, ...,
                                   selectionName = "Kolmogorov-Smirnov Test", verbose = 3)
  \S4method{KolmogorovSmirnovSelection}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
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
  \item{datasetName}{A name for the dataset used. Stored in the result.}                    
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
  Features are sorted in order of biggest distance to smallest. The top number of features
  is used in a classifier, to determine which number of features has the best resubstitution
  performance.
  
  Data tables which consist entirely of non-numeric data cannot be analysed. If \code{measurements}
  is an object of class \code{MultiAssayExperiment}, the factor of sample classes must be stored
  in the DataFrame accessible by the \code{colData} function with column name \code{"class"}.
}
\value{
  An object of class \code{\link{SelectResult}} or a list of such objects, if the classifier which
  was used for determining the specified performance metric made a number of prediction varieties.
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
  
  resubstituteParams <- ResubstituteParams(nFeatures = seq(10, 100, 10),
                                           performanceType = "balanced error",
                                           better = "lower")
  selected <- KolmogorovSmirnovSelection(genesMatrix, classes, "Example",
                             trainParams = TrainParams(naiveBayesKernel),
                             predictParams = PredictParams(NULL,
                                             getClasses = function(result) result),
                             resubstituteParams = resubstituteParams)
                             
  head(selected[["weighted=unweighted"]]@chosenFeatures)
  plotFeatureClasses(genesMatrix, classes, "Gene 13", dotBinWidth = 0.25,
                     xAxisLabel = expression(log[2](expression)))
  
}
