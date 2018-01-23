\name{DLDAinterface}
\alias{DLDAtrainInterface}
\alias{DLDAtrainInterface,matrix-method}
\alias{DLDAtrainInterface,DataFrame-method}
\alias{DLDAtrainInterface,MultiAssayExperiment-method}
\alias{DLDApredictInterface}
\alias{DLDApredictInterface,dlda,matrix-method}
\alias{DLDApredictInterface,dlda,DataFrame-method}
\alias{DLDApredictInterface,dlda,MultiAssayExperiment-method}
\title{
  An Interface for sparsediscrim Package's dlda Function
}
\description{
  \code{DLDAtrainInterface} generates a trained diagonal LDA classifier and \code{DLDApredictInterface}
  uses it to make predictions on a test dataset.
}
\usage{
  \S4method{DLDAtrainInterface}{matrix}(measurements, classes, ...)
  \S4method{DLDAtrainInterface}{DataFrame}(measurements, classes, ...)
  \S4method{DLDAtrainInterface}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
  \S4method{DLDApredictInterface}{matrix}(model, test, ...)
  \S4method{DLDApredictInterface}{DataFrame}(model, test, ...)
  \S4method{DLDApredictInterface}{MultiAssayExperiment}(model, test, targets = names(measurements), ...)  
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
                     \code{\link{MultiAssayExperiment}} containing the training data.
                     For a \code{matrix}, the rows are features, and the columns are samples.
                     If of type \code{\link{DataFrame}}, the dataset is subset to only those
                     features of type \code{integer}.}
  \item{classes}{Either a vector of class labels of class \code{\link{factor}} of the same length
                 as the number of samples in \code{measurements} or if the measurements are
                 of class \code{DataFrame} a character vector of length 1 containing the
                 column name in \code{measurement} is also permitted. Not used if \code{measurements}
                 is a \code{MultiAssayExperiment} object.}
  \item{test}{An object of the same class as \code{measurements} with no samples in common with
              \code{measurements} and the same number of features as it.}                 
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 integer variables from the clinical data table will be used.}
  \item{...}{Unused variables by the three top-level methods passed to the internal method
             which does the training or prediction (i.e. \code{verbose}).}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}  
}
\details{
  If \code{measurements} is an object of class \code{MultiAssayExperiment}, the factor of
  sample classes must be stored in the DataFrame accessible by the \code{colData} function
  with column name \code{"class"}.
}
\value{
  For \code{DLDAtrainInterface}, a trained DLDA classifier.
  For \code{DLDApredictInterface}, a result of type \code{list} with elements \code{class},
  \code{scores} and \code{posterior}, as created by \code{sparsediscrim}'s \code{predict} method.
}
\author{Dario Strbenac}
\examples{
  # Genes 76 to 100 have differential expression.
  genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
  genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
                                    c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
  classes <- factor(rep(c("Poor", "Good"), each = 25))
  colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix))
  rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
  selected <- rownames(genesMatrix)[91:100]
  trainingSamples <- c(1:20, 26:45)
  testingSamples <- c(21:25, 46:50)
  
  classifier <- DLDAtrainInterface(genesMatrix[selected, trainingSamples],
                                   classes[trainingSamples])
  DLDApredictInterface(classifier, genesMatrix[selected, testingSamples])[["class"]]
}
