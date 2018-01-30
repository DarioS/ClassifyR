\name{elasticNetGLMinterface}
\alias{elasticNetGLMinterface}
\alias{elasticNetGLMtrainInterface}
\alias{elasticNetGLMpredictInterface}
\alias{elasticNetGLMtrainInterface,matrix-method}
\alias{elasticNetGLMtrainInterface,DataFrame-method}
\alias{elasticNetGLMtrainInterface,MultiAssayExperiment-method}
\alias{elasticNetGLMpredictInterface,multnet,matrix-method}
\alias{elasticNetGLMpredictInterface,multnet,DataFrame-method}
\alias{elasticNetGLMpredictInterface,multnet,MultiAssayExperiment-method}
\title{
  An Interface for glmnet Package's glmnet Function
}
\description{
  An elastic net GLM classifier uses a penalty which is a combination of a lasso penality and
  a ridge penalty, scaled by a lambda value, to fit a sparse linear model to the data.
}
\usage{
  \S4method{elasticNetGLMtrainInterface}{matrix}(measurements, classes, ...)
  \S4method{elasticNetGLMtrainInterface}{DataFrame}(measurements, classes, ...)
  \S4method{elasticNetGLMtrainInterface}{MultiAssayExperiment}(measurements,
                                           targets = names(measurements), ...)
  \S4method{elasticNetGLMpredictInterface}{multnet,matrix}(model, test, ...)
  \S4method{elasticNetGLMpredictInterface}{multnet,DataFrame}(model, test, ...)
  \S4method{elasticNetGLMpredictInterface}{multnet,MultiAssayExperiment}(model, test,
                                           targets = names(measurements), ...)                                   
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
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method (e.g. \code{verbose}) or options
             that are used by the \code{link[glmnet]{glmnet}} function.}
}
\details{
  If \code{measurements} is an object of class \code{MultiAssayExperiment}, the factor of sample
  classes must be stored in the DataFrame accessible by the \code{colData} function with
  column name \code{"class"}.
  
  The value of the \code{"family"} parameter is set to \code{"multinomial"} so that
  multi-class classification is enabled.
}
\value{
  An object of type \code{randomForest}. The predictions of the test set samples are stored in
  the list element named \code{"predicted"} of the \code{"test"} element.
}
\author{Dario Strbenac}
\examples{
  if(require(glmnet))
  {
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
    
    trained <- elasticNetGLMtrainInterface(genesMatrix[, trainingSamples],
                                           classes[trainingSamples])
    elasticNetGLMpredictInterface(trained, genesMatrix[, testingSamples])
  }
}
