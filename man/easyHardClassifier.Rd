\name{easyHardClassifier}
\alias{easyHardClassifier}
\alias{easyHardClassifierTrain}
\alias{easyHardClassifierTrain,MultiAssayExperiment-method}
\alias{easyHardClassifierPredict}
\alias{easyHardClassifierPredict,EasyHardClassifier,MultiAssayExperiment-method}
\title{
  Two-stage Classification Using Easy-to-collect Data Set and Hard-to-collect data set.
}
\description{
  An alternative implementation to the previously published easy-hard classifier that doesn't do nested cross-validation for speed. In the first stage, each numeric variable is split on all possible midpoints between consecutive ordered values and the samples below the split and above the split are checked to see if they mostly belong to one class. Categorical varaibles are tabulated on factor levels and the count of samples in each class is determined. If any partitions of samples are pure for a class, based on a purity threshold, prediction rules are created. The samples not classified by any rule or classified to two or more classes the same number of times are left to be trained by the hard classifier.
}
\usage{
  \S4method{easyHardClassifierTrain}{MultiAssayExperiment}(measurements, easyDatasetID = "clinical", hardDatasetID = names(measurements)[1],
         featureSets = NULL, metaFeatures = NULL, minimumOverlapPercent = 80,
         easyClassifierParams = list(minCardinality = 10, minPurity = 0.9),
         hardClassifierParams = list(SelectParams(), TrainParams(), PredictParams()),
         characteristics = DataFrame(characteristic = c("Classifier Name", "Easy Dataset", "Hard Dataset"),
                                    value = c("Easy-Hard Classifier", easyDatasetID, hardDatasetID)), 
         verbose = 3)
  \S4method{easyHardClassifierPredict}{EasyHardClassifier,MultiAssayExperiment}(model, test, predictParams, verbose = 3)
}
\arguments{
  \item{measurements}{A \code{\link{MultiAssayExperiment}} object containing the data set The
                 sample classes must be in a column of the \code{DataFrame} accessed by \code{colData} named \code{"class"}}.
  \item{easyDatasetID}{The name of a data set in \code{measurements} or "clinical" to
  indicate the patient information in the column data be used.}
  \item{hardDatasetID}{The name of a data set in \code{measurements} different to the
  value of \code{easyDatasetID} to be used for classifying the samples not classified by the easy classifier.}
  \item{featureSets}{An object of type \code{\link{FeatureSetCollection}} which defines sets of
                     features or sets of edges.}               
  \item{metaFeatures}{Either \code{NULL} or a \code{DataFrame} which has meta-features
                      of the numeric data of interest.}
  \item{minimumOverlapPercent}{If \code{featureSets} stores sets of features, the minimum overlap
                               of feature IDs with \code{measurements} for a feature set to be retained
                               in the analysis. If \code{featureSets} stores sets of network edges,
                               the minimum percentage of edges with both vertex IDs found in \code{measurements}
                               that a set has to have to be retained in the analysis.}
  \item{easyClassifierParams}{A list of length 2 with names "minCardinality" and "minPurity". The first parameter specifies
  what the minimum number of samples after a split has to be and the second specifies the minimum proportion of samples in a partition belonging to a particular class.}
  \item{hardClassifierParams}{A list of objects defining the classification to do on the samples which were not predicted
  by the easy classifier. Objects must be of class \code{\link{TransformParams}}, \code{\link{SelectParams}},
  \code{\link{TrainParams}} or \code{\link{PredictParams}}.}
  \item{characteristics}{A \code{\link{DataFrame}} describing the characteristics
            of classification done. First column must be named \code{"charateristic"} and second column must
            be named \code{"value"}. If using wrapper functions for feature selection and classifiers
            in this package, the function names will automatically be
            generated and therefore it is not necessary to specify them.}
  \item{model}{A trained \code{\link{EasyHardClassifier}} object.}
  \item{test}{A \code{\link{MultiAssayExperiment}} object containing the test data.}
  \item{predictParams}{An object of class \code{\link{PredictParams}}. It specifies the classifier used to make
  the hard predictions.}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}             
}
\details{
  The easy classifier may be \code{NULL} if there are no rules that predicted the sample well using the easy data set. The hard classifier may be \code{NULL} if all of the samples could be predicted with rules generated using the easy data set or it will simply be a character if all or almost all of the remaining samples belong to one class.
}
\value{
  For \code{EasyHardClassifierTrain}, the trained two-stage classifier. For \code{EasyHardClassifierPredict},
  a factor vector of predicted classes.
}
\author{Dario Strbenac}
\references{
  Inspired by: Stepwise Classification of Cancer Samples Using Clinical and Molecular Data, Askar Obulkasim, Gerrit Meijer and Mark van de Wiel 2011,
  \emph{BMC Bioinformatics}, Volume 12 article 422, \url{https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-422}.
}
\examples{
  genesMatrix <- matrix(c(rnorm(90, 9, 1),
                          9.5, 9.4, 5.2, 5.3, 5.4, 9.4, 9.6, 9.9, 9.1, 9.8),
		        ncol = 10, byrow = TRUE)
  colnames(genesMatrix) <- paste("Sample", 1:10)
  rownames(genesMatrix) <- paste("Gene", 1:10)
  genders <- factor(c("Male", "Male", "Female", "Female", "Female",
                      "Female", "Female", "Female", "Female", "Female"))

  # Scenario: Male gender can predict the hard-to-classify Sample 1 and Sample 2.
  clinical <- DataFrame(age = c(31, 34, 32, 39, 33, 38, 34, 37, 35, 36),
                        gender = genders,
                        class = factor(rep(c("Poor", "Good"), each = 5)),
		        row.names = colnames(genesMatrix))
  dataset <- MultiAssayExperiment(ExperimentList(RNA = genesMatrix), clinical)
  selParams <- SelectParams(featureSelection = differentMeansSelection, selectionName = "Difference in Means",
                            resubstituteParams = ResubstituteParams(1:10, "balanced error", "lower"))
  trained <- easyHardClassifierTrain(dataset, easyClassifierParams = list(minCardinality = 2, minPurity = 0.9),
                                     hardClassifierParams = list(selParams, TrainParams(), PredictParams()))

  predictions <- easyHardClassifierPredict(trained, dataset, PredictParams())
}
