\name{runTests}
\alias{runTests}
\alias{runTestsEasyHard}
\alias{runTests,matrix-method}
\alias{runTests,DataFrame-method}
\alias{runTests,MultiAssayExperiment-method}
\alias{runTestsEasyHard,MultiAssayExperiment-method}
\title{Reproducibly Run Various Kinds of Cross-Validation}
\description{Enables doing classification schemes such as ordinary 10-fold,
             100 permutations 5-fold, and leave one out cross-validation.
             Processing in parallel is possible by leveraging the package \code{\link{BiocParallel}}.
             
             Pre-validation is possible and activated by specifying a list named \code{"prevaliated"} to
             \code{params}, which will use the functions specified in the list of parameters on the
             pre-validated data table. Other named items in the list correspond to other assays to be
             added as a pre-vaildated vector of the clinical table.
}
\usage{
  \S4method{runTests}{matrix}(measurements, classes, ...)
  \S4method{runTests}{DataFrame}(measurements, classes, featureSets = NULL, metaFeatures = NULL,
         minimumOverlapPercent = 80, datasetName, classificationName,
         validation = c("permute", "leaveOut", "fold"),
         permutePartition = c("fold", "split"),
         permutations = 100, percent = 25, folds = 5, leave = 2,
         seed, parallelParams = bpparam(),
            params = list(SelectParams(), TrainParams(), PredictParams()), verbose = 1)
  \S4method{runTests}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
  \S4method{runTestsEasyHard}{MultiAssayExperiment}(measurements, easyDatasetID = "clinical", hardDatasetID = names(measurements)[1],
                   featureSets = NULL, metaFeatures = NULL, minimumOverlapPercent = 80,
                   datasetName = NULL, classificationName = "Easy-Hard Classifier", 
                   validation = c("permute", "leaveOut", "fold"),
                   permutePartition = c("fold", "split"),
                   permutations = 100, percent = 25, folds = 5, leave = 2,
                   seed, parallelParams = bpparam(), ..., verbose = 1)  
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
                      \code{\link{MultiAssayExperiment}} containing the training data.
                      For a \code{matrix}, the rows are features, and the columns are samples.
                      The sample identifiers must be present as column names of the \code{matrix}
                      or the row names of the \code{DataFrame}. If pre-validation is activated by naming
                      one of the lists in \code{params} \code{"prevalidated"}, this variable must be of type
                      \code{\link{MultiAssayExperiment}}.}
  \item{classes}{Either a vector of class labels of class \code{\link{factor}} of the same length
                 as the number of samples in \code{measurements} or if the measurements are
                 of class \code{DataFrame} a character vector of length 1 containing the
                 column name in \code{measurement} is also permitted. Not used if \code{measurements}
                 is a \code{MultiAssayExperiment} object.}
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
  \item{...}{For \code{runTests}, variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method. For \code{runTestsEasyHard}, \code{easyClassifierParams}
             and \code{hardClassifierParams} to be passed to \code{\link{easyHardClassifierTrain}}}        
  \item{datasetName}{A name associated with the data set used.}
  \item{classificationName}{A name associated with the classification.}
  \item{validation}{Default: \code{"permute"}. \code{"permute"} for repeated permuting.
                    \code{"leaveOut"} for leaving all possible combinations of k samples as test samples.
                    \code{"fold"} for folding of the data set (no resampling).}
  \item{permutePartition}{Default: \code{"fold"}. Either \code{"fold"} or \code{"split"}. Only applicable if
                          \code{validation} is \code{"permute"}. If \code{"fold"}, then the samples are
                          split into folds and in each iteration one is used as the test set. If
                          \code{"split"}, the samples are split into two groups, the sizes being based on
                          the \code{percent} value. One group is used as the training set,
                          the other is the test set.
                         }                    
  \item{permutations}{Default: 100. Relevant when permuting is used. The number of times to do reordering
                      of the samples before splitting or folding them.}
  \item{percent}{Default: 25. Used when permutation with the split method is chosen. The percentage of
                 samples to be in the test set.}        
  \item{folds}{Default: 5. Relevant when repeated permutations are done and \code{permutePartition} is
               set to \code{"fold"} or when \code{validation} is set to \code{"fold"}. The number of folds
               to break the data set into. Each fold is used once as the test set.}
  \item{leave}{Default: 2. Relevant when leave-k-out cross-validation is used. The number of samples to
               leave for testing.}
  \item{seed}{The random number generator used for repeated resampling will use this seed,
              if it is provided. Allows reproducibility of repeated usage on the same input data.}
  \item{parallelParams}{An object of class \code{\link{MulticoreParam}} or \code{\link{SnowParam}}.}          
  \item{params}{A \code{list} of objects of class of \code{\link{TransformParams}},
               \code{\link{SelectParams}}, \code{\link{TrainParams}} or \code{\link{PredictParams}}.
               The order they are in the list determines the order in which the stages
               of classification are done in. It may also be a list of such lists for pre-validation.
               In that case, each list must be named and one of them must be named \code{"prevalidated"},
               which specifies the functions to use on the pre-validated data table.
               }
  \item{easyDatasetID}{The name of a data set in \code{measurements} or "clinical" to
  indicate the patient information in the column data be used.}
  \item{hardDatasetID}{The name of a data set in \code{measurements} different to the
  value of \code{easyDatasetID} to be used for classifying the samples not classified by the easy classifier.}               
  \item{verbose}{Default: 1. A number between 0 and 3 for the amount of progress messages to give.
                 A higher number will produce more messages as more lower-level functions
                 print messages.}
}
\value{
  If the predictor function made a single prediction per sample, then an object of class
  \code{\link{ClassifyResult}}. If the predictor function made a set of predictions, then
  a list of such objects.
}
\author{Dario Strbenac}

\examples{
  #if(require(sparsediscrim))
  #{
    data(asthma)
    
    resubstituteParams <- ResubstituteParams(nFeatures = seq(5, 25, 5),
                                         performanceType = "balanced error",
                                         better = "lower")
    runTests(measurements, classes, datasetName = "Asthma",
             classificationName = "Different Means", permutations = 5,
             params = list(SelectParams(differentMeansSelection, "t Statistic",
                                        resubstituteParams = resubstituteParams),
                           TrainParams(DLDAtrainInterface),
                           PredictParams(DLDApredictInterface)
                           )
             )
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
  easyHardCV <- runTestsEasyHard(dataset, datasetName = "Test Data", classificationName = "Easy-Hard",
                                 easyClassifierParams = list(minCardinality = 2, minPurity = 0.9),
                                 hardClassifierParams = list(selParams, TrainParams(), PredictParams()),
                                 validation = "leaveOut", leave = 1)  
}
