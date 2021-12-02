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
  \S4method{runTests}{DataFrame}(measurements, classes, crossValParams = CrossValParams(),
                   modellingParams = ModellingParams(), characteristics = DataFrame(), verbose = 1)
  \S4method{runTests}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
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
  \item{crossValParams}{An object of class \code{\link{CrossValParams}}, specifying the kind of
  cross-validation to be done.}
  \item{modellingParams}{An object of class \code{\link{ModellingParams}}, specifying the class rebalancing, transformation (if any),
  feature selection (if any), training and prediction to be done on the data set.}  
  \item{characteristics}{A \code{\link{DataFrame}} describing the characteristics
            of the classification used. First column must be named \code{"charateristic"} and second column
            must be named \code{"value"}. Useful for automated plot annotation by plotting functions within this package.
            Transformation, selection and prediction functions provided by this package will cause the characteristics
            to be automatically determined and this can be left blank.}                               
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 numeric variables from the clinical data table will be used.}
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method.}
  \item{verbose}{Default: 1. A number between 0 and 3 for the amount of progress messages to give.
                 A higher number will produce more messages as more lower-level functions
                 print messages.}
}
\value{
  An object of class \code{\link{ClassifyResult}}.
}
\author{Dario Strbenac}

\examples{
  #if(require(sparsediscrim))
  #{
    data(asthma)
    
    CVparams <- CrossValParams(permutations = 5)
    tuneList <- list(nFeatures = seq(5, 25, 5), performanceType = "Balanced Error")
    selectParams <- SelectParams(differentMeansRanking, tuneParams = tuneList)
    modellingParams <- ModellingParams(selectParams = selectParams)
    runTests(measurements, classes, CVparams, modellingParams,
             DataFrame(characteristic = c("Dataset Name", "Classifier Name"),
                       value = c("Asthma", "Different Means"))
             )
  #}
}
