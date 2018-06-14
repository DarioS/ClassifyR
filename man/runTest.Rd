\name{runTest}
\alias{runTest}
\alias{runTest,matrix-method}
\alias{runTest,DataFrame-method}
\alias{runTest,MultiAssayExperiment-method}
\title{Perform a Single Classification}
\description{For a dataset of features and samples, the classification process is run.
             It consists of data transformation, feature selection, classifier training
             and testing (prediction of samples not used in training).
       }
\usage{
  \S4method{runTest}{matrix}(measurements, classes, ...)
  \S4method{runTest}{DataFrame}(measurements, classes, datasetName, classificationName,
            training, testing, params = list(SelectParams(), TrainParams(), PredictParams()),
        verbose = 1, .iteration = NULL)
  \S4method{runTest}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
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
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 numeric variables from the clinical data table will be used.}
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method.}
  \item{datasetName}{A name associated with the dataset used.}
  \item{classificationName}{A name associated with the classification.}  
  \item{training}{A vector which specifies the training samples.}
  \item{testing}{A vector which specifies the test samples.}    
  \item{params}{A \code{\link{list}} of objects of class of \code{\link{TransformParams}},
               \code{\link{SelectParams}}, \code{\link{TrainParams}}, or \code{\link{PredictParams}}.
               The order they are in the list determines the order in which the stages
               of classification are done in.}
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
  if(require(sparsediscrim))
  {
    data(asthma)
    resubstituteParams <- ResubstituteParams(nFeatures = seq(5, 25, 5),
                                         performanceType = "balanced error",
                                         better = "lower")
    runTest(measurements, classes, "Asthma", "Different Means",
            params = list(SelectParams(limmaSelection, "Moderated t Statistic",
                                       resubstituteParams = resubstituteParams),
                          TrainParams(DLDAtrainInterface),
                          PredictParams(DLDApredictInterface,
                                        getClasses = function(result) result[["class"]])),
            training = (1:ncol(measurements)) \%\% 2 == 0,
            testing = (1:ncol(measurements)) \%\% 2 != 0)
  }
}
