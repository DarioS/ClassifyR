\name{networkCorrelationsSelection}
\alias{networkCorrelationsSelection}
\alias{networkCorrelationsSelection,matrix-method}
\alias{networkCorrelationsSelection,DataFrame-method}
\alias{networkCorrelationsSelection,MultiAssayExperiment-method}
\title{Selection of Differentially Correlated Hub Sub-networks}
\description{Ranks sub-networks by largest within-class to between-class
             correlation variability and chooses the sub-networks which have
             the best resubstitution performance.
}
\usage{
  \S4method{networkCorrelationsSelection}{matrix}(measurements, classes, metaFeatures = NULL, ...)
  \S4method{networkCorrelationsSelection}{DataFrame}(measurements, classes, metaFeatures = NULL,
                  featureSets, datasetName, trainParams, predictParams, resubstituteParams,
                  selectionName = "Differential Correlation of Sub-networks", verbose = 3)
  \S4method{networkCorrelationsSelection}{MultiAssayExperiment}(measurements, target = NULL, metaFeatures = NULL, ...)
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
  \item{metaFeatures}{A \code{DataFrame} with the same number of samples as the numeric table
                      of interest. The number of derived features in this table will be different
                      to the original input data table. The command \code{mcols(metaFeatures)}
                      must return a \code{DataFrame} which has an "original" column with as many rows
                      as there are meta-features and specifies the feature which the meta-feature is
                      originally derived from (e.g. network name).}
  \item{featureSets}{A object of type \code{\link{FeatureSetCollection}}. The \code{sets} slot must contain
                     a list of two-column matrices with each row corresponding to a binary interaction. Such
                     sub-networks may be determined by a community detection algorithm. This will be used
                     to determine which features belong to which sub-networks before calculating a
                     statistic for each sub-network.}
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
  The selection of sub-networks is based on the average difference in correlations between
  each pair of interactors, considering the samples within each class separately. Such differences
  of correlations within each of the two classes are scaled by the average difference of correlations
  within each class.
  
  More formally, let \eqn{C_{i,j}}{C_{i,j}} be the correlation of the
  \eqn{j^{\textnormal{th}}}{j-th} edge using all samples belonging to to class \eqn{i}.
  Then, let \eqn{\overline{C}_{i\bullet}}{mean(C_{i,*})} be defined as
  \eqn{\frac{C_{i,1} + C_{i,2} + ... + C_{i,e}}{n}}{(C_{i,1} + C_{i,2} + ... + C_{i,e}) / e}
  where \eqn{e}{e} is the number of edges in the sub-network being considered. Also, let
  \eqn{\overline{C}_{\bullet\bullet}}{mean(C_{*,*})}, the average overall correlation, be
  \eqn{\frac{C_{1\bullet} + C_{2\bullet}}{2}}{(C_{1,*} + C_{1,*}) / 2}. Then, the
  between-class sum-of-squares (BSS) is \eqn{\sum_{i=1}^{2}e(\overline{C}_{i\bullet} - \overline{C}_{\bullet\bullet})^2}{(mean(C_{1,*}) - mean(C_{*,*})^2 + (mean(C_{2,*}) - mean(C_{*,*})^2}. Also the within-class
  sum-of-squares (WSS) is \eqn{\sum_{i=1}^{2}\sum_{j=1}^{e}(C_{i,j} - \overline{C}_{i\bullet})^2}{sum(sum((C_{i,j} - mean(C_{i,*}))^2, j is 1 to e), i is 1 to 2)}. The sub-networks are ranked in decreasing 
  order of \eqn{\frac{BSS}{WSS}}{BSS/WSS}.
  
  The classifier specified by \code{trainParams} and \code{predictParams} is used to calculate resubtitution
  error rates using the transformation of the data set provided by \code{metaFeatures}.
  The set of top-ranked sub-networks which give the lowest resubstitution error rate are finally selected.
  
  Data tables which consist entirely of non-numeric data cannot be analysed. If \code{measurements}
  is an object of class \code{MultiAssayExperiment}, the factor of sample classes must be stored
  in the DataFrame accessible by the \code{colData} function with column name \code{"class"}.
}
\value{
  An object of class \code{\link{SelectResult}} or a list of such objects, if the classifier which was
  used for determining the specified performance metric made a number of prediction varieties.
}
\author{Dario Strbenac}
\seealso{
\code{\link{interactorDifferences}} for an example of a function which can turn the measurements into
  meta-features for classification.
}
\references{
  Network-based biomarkers enhance classical approaches to prognostic gene expression signatures,
  Rebecca L Barter, Sarah-Jane Schramm, Graham J Mann and Yee Hwa Yang, 2014,
  \emph{BMC Systems Biology}, Volume 8 Supplement 4 Article S5, \url{https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-8-S4-S5}.
}
\examples{
  networksList <- list(`A Hub` = matrix(c('A', 'A', 'A', 'B', 'C', 'D'), ncol = 2),
                       `G Hub` = matrix(c('G', 'G', 'G', 'H', 'I', 'J'), ncol = 2))
  netSets <- FeatureSetCollection(networksList)
                           
  # Differential correlation for sub-network with hub A.                                           
  measurements <- matrix(c(5.7, 10.1, 6.9, 7.7, 8.8, 9.1, 11.2, 6.4, 7.0, 5.5,
                           5.6, 9.6, 7.0, 8.4, 10.8, 12.2, 8.1, 5.7, 5.4, 12.1,
                           4.5, 9.0, 6.9, 7.0, 7.3, 6.9, 7.8, 7.9, 5.7, 8.7,
                           8.1, 10.6, 7.4, 7.1, 10.4, 6.1, 7.3, 2.7, 11.0, 9.1,
                           round(rnorm(60, 8, 1), 1)), ncol = 10, byrow = TRUE)
  classes <- factor(rep(c("Good", "Poor"), each = 5))
                         
  rownames(measurements) <- LETTERS[1:10]
  colnames(measurements) <- names(classes) <- paste("Patient", 1:10)
  
  Idifferences <- interactorDifferences(measurements, netSets)

  # The features are sub-networks and there are only two in this example.
  resubstituteParams <- ResubstituteParams(nFeatures = 1:2,
                                performanceType = "balanced error")
  
  predictParams <- PredictParams(NULL)
  networkCorrelationsSelection(measurements, classes, metaFeatures = Idifferences,
                               featureSets = netSets,
                               trainParams = TrainParams(naiveBayesKernel),
                               predictParams = predictParams,
                               resubstituteParams = resubstituteParams)
}
