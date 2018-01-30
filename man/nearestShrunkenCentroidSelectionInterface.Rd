\name{nearestShrunkenCentroidSelectionInterface}
\alias{nearestShrunkenCentroidSelectionInterface}
\alias{nearestShrunkenCentroidSelectionInterface,matrix-method}
\alias{nearestShrunkenCentroidSelectionInterface,DataFrame-method}
\alias{nearestShrunkenCentroidSelectionInterface,MultiAssayExperiment-method}

\title{Interface for \code{pamr.listgenes} Function from \code{pamr} CRAN Package}
\description{Restructures variables from ClassifyR framework to be compatible with \code{\link[pamr]{pamr.listgenes}} definition.}
\usage{
  \S4method{nearestShrunkenCentroidSelectionInterface}{matrix}(measurements, classes, ...)
  \S4method{nearestShrunkenCentroidSelectionInterface}{DataFrame}(measurements, classes, datasetName,
                   trained, ..., selectionName = "Shrunken Centroids", verbose = 3)
  \S4method{nearestShrunkenCentroidSelectionInterface}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
                      \code{\link{MultiAssayExperiment}} containing the training data.
                      For a \code{matrix}, the rows are features, and the columns are samples.}
  \item{classes}{Either a vector of class labels of class \code{\link{factor}} of the same length
                 as the number of samples in \code{measurements} or if the measurements are
                 of class \code{DataFrame} a character vector of length 1 containing the
                 column name in \code{measurement} is also permitted. Not used if \code{measurements}
                 is a \code{MultiAssayExperiment} object.}
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 numeric variables from the clinical data table will be used.}                 
  \item{datasetName}{A name for the dataset used. Stored in the result.}
  \item{trained}{The output of \code{\link{nearestShrunkenCentroidTrainInterface}},
                      which is identical to the output of \code{\link[pamr]{pamr.listgenes}}.}
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method
             which are passed into and used by the \code{DataFrame} method or extra arguments
             passed to \code{\link[pamr]{pamr.listgenes}}.}
  \item{selectionName}{A name to identify this selection method by. Stored in the result.}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  This function is an interface between the ClassifyR framework and \code{\link[pamr]{pamr.listgenes}}.
  
  The set of features chosen is the obtained by considering the range of thresholds provided to \code{\link{nearestShrunkenCentroidTrainInterface}} and using the threshold that obtains the lowest cross-validation error rate on the training set.
}
\value{
  An object of class \code{\link{SelectResult}}. The \code{rankedFeatures} slot will be empty.
}
\author{Dario Strbenac}
\seealso{\code{\link[pamr]{pamr.listgenes}} for the function that was interfaced to.}
\examples{
  if(require(pamr))
  {
    # Genes 76 to 100 have differential expression.
    genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(100, 9, 1)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn)
                                 c(rnorm(75, 9, 1), rnorm(25, 14, 1))))
    rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))                                 
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    
    trained <- nearestShrunkenCentroidTrainInterface(genesMatrix, classes)
    selected <- nearestShrunkenCentroidSelectionInterface(genesMatrix, classes, "Example",
                                                          trained)
                                                          
    selected@chosenFeatures                                                         
  }
}
