\name{runTest}
\alias{runTest}
\alias{runTest,matrix-method}
\alias{runTest,DataFrame-method}
\alias{runTest,MultiAssayExperiment-method}
\title{Perform a Single Classification}
\description{For a data set of features and samples, the classification process is run.
             It consists of data transformation, feature selection, classifier training
             and testing.
       }
\usage{
  \S4method{runTest}{matrix}(measurements, classes, ...)
  \S4method{runTest}{DataFrame}(measurements, classes, training, testing, crossValParams = CrossValParams(),
         modellingParams = ModellingParams(), characteristics = DataFrame(), verbose = 1, .iteration = NULL)
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
  \item{crossValParams}{An object of class \code{\link{CrossValParams}}, specifying the kind of
  cross-validation to be done, if nested cross-validation is used to tune any parameters.}
  \item{modellingParams}{An object of class \code{\link{ModellingParams}}, specifying the class rebalancing, transformation (if any),
  feature selection (if any), training and prediction to be done on the data set.}               
                               
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 numeric variables from the clinical data table will be used.}
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method.}
  \item{characteristics}{A \code{\link{DataFrame}} describing the characteristics
            of the classification used. First column must be named \code{"charateristic"} and second column
            must be named \code{"value"}. Useful for automated plot annotation by plotting functions within this package.
            Transformation, selection and prediction functions provided by this package will cause the characteristics
            to be automatically determined and this can be left blank.}             
  \item{training}{A vector which specifies the training samples.}
  \item{testing}{A vector which specifies the test samples.}    
  \item{verbose}{Default: 1. A number between 0 and 3 for the amount of progress messages to give.
                 A higher number will produce more messages as more lower-level functions
                 print messages.}
  \item{.iteration}{Not to be set by a user. This value is used to keep track of the cross-validation
                   iteration, if called by \code{\link{runTests}}.}                
}
\details{
  This function only performs one classification and prediction. See \code{\link{runTests}}
  for a driver function that enables a number of different cross-validation schemes to be applied
  and uses this function to perform each iteration.
}
\value{
  If called directly by the user rather than being used internally by \code{\link{runTests}}, a
  \code{\link{ClassifyResult}} object. Otherwise a list of different aspects of the result which is
  passed back to \code{\link{runTests}}.
}
\author{Dario Strbenac}

\examples{
  #if(require(sparsediscrim))
  #{
    data(asthma)
    tuneList <- list(nFeatures = seq(5, 25, 5), performanceType = "Balanced Error")
    selectParams <- SelectParams(limmaRanking, tuneParams = tuneList)
    modellingParams <- ModellingParams(selectParams = selectParams)
    runTest(measurements, classes, modellingParams = modellingParams,
            training = seq(1, ncol(measurements), 2), testing = seq(2, ncol(measurements), 2))
  #}
}
