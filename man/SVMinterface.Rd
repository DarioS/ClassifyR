\name{SVMinterface}
\alias{SVMtrainInterface}
\alias{SVMtrainInterface,matrix-method}
\alias{SVMtrainInterface,DataFrame-method}
\alias{SVMtrainInterface,MultiAssayExperiment-method}
\alias{SVMpredictInterface}
\alias{SVMpredictInterface,svm,matrix-method}
\alias{SVMpredictInterface,svm,DataFrame-method}
\alias{SVMpredictInterface,svm,MultiAssayExperiment-method}
\title{
  An Interface for e1071 Package's Support Vector Machine Classifier.
}
\description{
  \code{SVMtrainInterface} generates a trained SVM classifier and \code{SVMpredictInterface}
  uses it to make predictions on a test data set.
}
\usage{
  \S4method{SVMtrainInterface}{matrix}(measurements, classes, ...)
  \S4method{SVMtrainInterface}{DataFrame}(measurements, classes, ..., verbose = 3)
  \S4method{SVMtrainInterface}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
  \S4method{SVMpredictInterface}{svm,matrix}(model, test, ...)
  \S4method{SVMpredictInterface}{svm,DataFrame}(model, test, classes = NULL, 
                                 returnType = c("both", "class", "score"),
                                 verbose = 3)
  \S4method{SVMpredictInterface}{svm,MultiAssayExperiment}(model, test, targets = names(test), ...)  
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
                     \code{\link{MultiAssayExperiment}} containing the training data.
                     For a \code{matrix}, the rows are features, and the columns are samples.
                     If of type \code{\link{DataFrame}}, the data set is subset to only those
                     features of type \code{integer}.}
  \item{classes}{Either a vector of class labels of class \code{\link{factor}} of the same length
                 as the number of samples in \code{measurements} or if the measurements are
                 of class \code{DataFrame} a character vector of length 1 containing the
                 column name in \code{measurement} is also permitted. Not used if \code{measurements}
                 is a \code{MultiAssayExperiment} object.}
  \item{returnType}{Default: \code{"both"}. Either \code{"class"}, \code{"score"} or \code{"both"}. Sets the return value
                    from the prediction to either a vector of class labels, score for a sample belonging
                    to the second class, as determined by the factor levels, or both labels and scores
                    in a \code{data.frame}.}                 
  \item{test}{An object of the same class as \code{measurements} with no samples in common with
              \code{measurements} and the same number of features as it.  Also, if a \code{DataFrame},
              the \code{class} column must be absent.}                 
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 numeric variables from the clinical data table will be used.}
  \item{model}{A fitted model as returned by \code{SVMtrainInterface}.}                 
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method (e.g. \code{verbose}) or options
             that are used by the \code{\link[e1071]{svm}} function.}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}  
}
\details{
  If \code{measurements} is an object of class \code{MultiAssayExperiment}, the factor of
  sample classes must be stored in the DataFrame accessible by the \code{colData} function
  with column name \code{"class"}.
}
\value{
  For \code{SVMtrainInterface}, a trained SVM classifier of type \code{svm}.
  For \code{SVMpredictInterface}, either a factor vector of predicted classes, a vector of secores for the
  second class, or a table of both the class labels and second class scores, depending on the setting
  of \code{returnType}.
}
\author{Dario Strbenac}
\examples{
  if(require(e1071))
  {
    # Genes 76 to 100 have differential expression.
    genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
                                      c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix))
    rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
    trainingSamples <- c(1:20, 26:45)
    testingSamples <- c(21:25, 46:50)
    
    classifier <- SVMtrainInterface(genesMatrix[, trainingSamples],
                                     classes[trainingSamples], kernel = "linear")
    SVMpredictInterface(classifier, genesMatrix[, testingSamples])
  }
}
