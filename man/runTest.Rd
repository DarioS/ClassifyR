\name{runTest}
\alias{runTest}
\alias{runTestEasyHard}
\alias{runTest,matrix-method}
\alias{runTest,DataFrame-method}
\alias{runTest,MultiAssayExperiment-method}
\alias{runTestEasyHard,MultiAssayExperiment-method}
\title{Perform a Single Classification}
\description{For a data set of features and samples, the classification process is run.
             It consists of data transformation, feature selection, classifier training
             and testing.
       }
\usage{
  \S4method{runTest}{matrix}(measurements, classes, ...)
  \S4method{runTest}{DataFrame}(measurements, classes,
            balancing = c("downsample", "upsample", "none"), featureSets = NULL, metaFeatures = NULL,
            minimumOverlapPercent = 80, datasetName, classificationName,
            training, testing, params = list(SelectParams(), TrainParams(), PredictParams()),
        verbose = 1, .iteration = NULL)
  \S4method{runTest}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
  \S4method{runTestEasyHard}{MultiAssayExperiment}(measurements, balancing = c("downsample", "upsample", "none"),
                   easyDatasetID = "clinical", hardDatasetID = names(measurements)[1],
                   featureSets = NULL, metaFeatures = NULL, minimumOverlapPercent = 80,
                   datasetName = NULL, classificationName = "Easy-Hard Classifier", training, testing, ..., verbose = 1, .iteration = NULL)
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
                      \code{\link{MultiAssayExperiment}} containing the training data.
                      For a \code{matrix}, the rows are features, and the columns are samples.
                      The sample identifiers must be present as column names of the \code{matrix}
                      or the row names of the \code{DataFrame}.}
  \item{classes}{Either a vector of class labels of class \code{\link{factor}} of the same length
                 as the number of samples in \code{measurements} or if the measurements are
                 of class \code{DataFrame} a character vector of length 1 containing the
                 column name in \code{measurement} is also permitted. Not used if \code{measurements}
                 is a \code{MultiAssayExperiment} object.}
  \item{balancing}{A character vector of length 1 specifying how to balance samples in the training set
                   so that each class ends up with an equal number of samples. \code{"downsample"}
                   samples from each class with more samples than the smallest class so that each of
                   those classes ends up with the same number of samples as the smallest class.
                   \code{"upsample"} causes sampling with replacement to be done for every class that
                   has less samples than the largest class to make all classes have equal size. Lastly,
                   \code{"none"} means that no downsampling or upsampling is done and the
                   class imbalance is preserved.}                 
  \item{featureSets}{An object of type \code{\link{FeatureSetCollection}} which defines sets of
                     features or sets of edges.}                 
  \item{metaFeatures}{Either \code{NULL} or a \code{DataFrame} which has meta-features
                      of the numeric data of interest.}
  \item{minimumOverlapPercent}{If \code{featureSets} stores sets of features, the minimum overlap
                               of feature IDs with \code{measurements} for a feature set to be retained
                               in the analysis. If \code{featureSets} stores sets of network edges,
                               the minimum percentage of edges with both vertex IDs found in \code{measurements}
                               that a set has to have to be retained in the analysis.}                      
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 numeric variables from the clinical data table will be used.}
  \item{...}{For \code{runTest}, variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method. For \code{runTestEasyHard}, \code{easyClassifierParams}
             and \code{hardClassifierParams} to be passed to \code{\link{easyHardClassifierTrain}}.}
  \item{datasetName}{A name associated with the data set used.}
  \item{classificationName}{A name associated with the classification.}  
  \item{training}{A vector which specifies the training samples.}
  \item{testing}{A vector which specifies the test samples.}    
  \item{params}{A \code{\link{list}} of objects of class of \code{\link{TransformParams}},
               \code{\link{SelectParams}}, \code{\link{TrainParams}}, or \code{\link{PredictParams}}.
               The order they are in the list determines the order in which the stages
               of classification are done in.}
  \item{easyDatasetID}{The name of a data set in \code{measurements} or "clinical" to
  indicate the patient information in the column data be used.}
  \item{hardDatasetID}{The name of a data set in \code{measurements} different to the
  value of \code{easyDatasetID} to be used for classifying the samples not classified by the easy classifier.}
  \item{verbose}{Default: 1. A number between 0 and 3 for the amount of progress messages to give.
                 A higher number will produce more messages as more lower-level functions
                 print messages.}
  \item{.iteration}{Not to be set by a user. This value is used to keep track of the cross-validation
                   iteration, if called by \code{\link{runTests}}.}                
}
\details{
  This function only performs one classification and prediction. See \code{\link{runTests}}
  for a driver function that enables a number of different cross-validation schemes to be applied
  and uses this function to perform each iteration. \code{datasetName} and \code{classificationName}
  need to be provided.
}
\value{
  If called directly by the user rather than being used internally by \code{runTests}, a
  \code{SelectResult} object.
}
\author{Dario Strbenac}

\examples{
  #if(require(sparsediscrim))
  #{
    data(asthma)
    resubstituteParams <- ResubstituteParams(nFeatures = seq(5, 25, 5),
                                         performanceType = "balanced error",
                                         better = "lower")
    runTest(measurements, classes,
            datasetName = "Asthma", classificationName = "Different Means",
            params = list(SelectParams(limmaSelection, "Moderated t Statistic",
                                       resubstituteParams = resubstituteParams),
                          TrainParams(DLDAtrainInterface),
                          PredictParams(DLDApredictInterface)
                          ),
            training = (1:ncol(measurements)) \%\% 2 == 0,
            testing = (1:ncol(measurements)) \%\% 2 != 0)
  #}
  
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
  easyHardCV <- runTestEasyHard(dataset, datasetName = "Test Data", classificationName = "Easy-Hard", training = paste("Sample", 1:10), testing = paste("Sample", 1:10),
                                easyClassifierParams = list(minCardinality = 2, minPurity = 0.9),
                                hardClassifierParams = list(selParams, TrainParams(), PredictParams())
                               )
}
