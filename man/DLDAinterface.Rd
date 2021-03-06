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
  uses it to make predictions on a test data set.
}
\usage{
  \S4method{DLDAtrainInterface}{matrix}(measurements, classes, ...)
  \S4method{DLDAtrainInterface}{DataFrame}(measurements, classes, verbose = 3)
  \S4method{DLDAtrainInterface}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
  \S4method{DLDApredictInterface}{dlda,matrix}(model, test, ...)
  \S4method{DLDApredictInterface}{dlda,DataFrame}(model, test, returnType = c("class", "score", "both"),
                                  verbose = 3)
  \S4method{DLDApredictInterface}{dlda,MultiAssayExperiment}(model, test, targets = names(test), ...)  
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
  \item{model}{A fitted model as returned by \code{DLDAtrainInterface}.}
  \item{test}{An object of the same class as \code{measurements} with no samples in common with
              \code{measurements} and the same number of features as it. Also, if a \code{DataFrame},
              the \code{class} column must be absent.}                 
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 integer variables from the clinical data table will be used.}
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method (e.g. \code{verbose}).}
  \item{returnType}{Default: \code{"class"}. Either \code{"class"}, \code{"score"} or \code{"both"}.
                    Sets the return value from the prediction to either a vector of class labels,
                    matrix of scores for each class, or both labels and scores in a \code{data.frame}.}   
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
  For \code{DLDApredictInterface}, either a factor vector of predicted classes, a matrix of
  scores for each class, or a table of both the class labels and class scores, depending on
  the setting of \code{returnType}.
}
\author{Dario Strbenac}
\examples{
  # if(require(sparsediscrim)) Package currently removed from CRAN.
  #{
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
    DLDApredictInterface(classifier, genesMatrix[selected, testingSamples])
  #}
}
