\name{runTests}
\alias{runTests}
\alias{runTests,matrix-method}
\alias{runTests,DataFrame-method}
\alias{runTests,MultiAssayExperiment-method}
\title{Reproducibly Run Various Kinds of Cross-Validation}
\description{Enables doing classification schemes such as ordinary 10-fold,
             100 permutations 5-fold, and leave one out cross-validation.
             Processing in parallel is possible by leveraging the package \code{\link{BiocParallel}}.
}
\usage{
  \S4method{runTests}{matrix}(measurements, classes, ...)
  \S4method{runTests}{DataFrame}(measurements, classes, metaFeatures = NULL, datasetName, classificationName,
         validation = c("permute", "leaveOut", "fold"),
         permutePartition = c("fold", "split"),
         permutations = 100, percent = 25, folds = 5, leave = 2,
         seed, parallelParams = bpparam(),
            params = list(SelectParams(), TrainParams(), PredictParams()), verbose = 1)
  \S4method{runTests}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
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
  \item{metaFeatures}{Either \code{NULL} or a \code{DataFrame} which has meta-features
                      of the numeric data of interest.}                 
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 numeric variables from the clinical data table will be used.}
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method.}        
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
               of classification are done in.}
  \item{verbose}{Default: 1. A number between 0 and 3 for the amount of progress messages to give.
                 A higher number will produce more messages as more lower-level functions
                 print messages.}
}
\value{
  If the predictor function made a single prediction, then an object of class
  \code{\link{ClassifyResult}}. If the predictor function made a set of predictions, then
  a list of such objects.
}
\author{Dario Strbenac}

\examples{
  if(require(sparsediscrim))
  {
    data(asthma)
    
    resubstituteParams <- ResubstituteParams(nFeatures = seq(5, 25, 5),
                                         performanceType = "balanced error",
                                         better = "lower")
    runTests(measurements, classes, "Asthma", "Different Means",
             permutations = 5,
             params = list(SelectParams(tTestSelection, "t Statistic",
                                        resubstituteParams = resubstituteParams),
                           TrainParams(DLDAtrainInterface),
                           PredictParams(DLDApredictInterface,
                                         getClasses = function(result) result[["class"]])))
  }
}
