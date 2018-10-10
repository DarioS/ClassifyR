\name{classifyInterface}
\alias{classifyInterface}
\alias{classifyInterface,matrix-method}
\alias{classifyInterface,DataFrame-method}
\alias{classifyInterface,MultiAssayExperiment-method}
\title{
  An Interface for PoiClaClu Package's Classify Function
}
\description{
  More details of Poisson LDA are available in the documentation of \code{\link[PoiClaClu]{Classify}}.
}
\usage{
  \S4method{classifyInterface}{matrix}(measurements, classes, test, ...)
  \S4method{classifyInterface}{DataFrame}(measurements, classes, test, ..., returnType = c("class", "score", "both"), verbose = 3)
  \S4method{classifyInterface}{MultiAssayExperiment}(measurements, test, targets = names(measurements), ...)
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
  \item{test}{An object of the same class as \code{measurements} with no samples in common with
              \code{measurements} and the same number of features as it.}                 
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 integer variables from the clinical data table will be used.}
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method or parameters that
             \code{\link[PoiClaClu]{Classify}} can accept.}
  \item{returnType}{Default: \code{"class"}. Either \code{"class"}, \code{"score"} or \code{"both"}.
                    Sets the return value from the prediction to either a vector of class labels,
                    matrix of scores for each class, or both labels and scores in a \code{data.frame}.}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}    
}
\details{
  Data tables which consist entirely of non-integer data cannot be analysed. If \code{measurements}
  is an object of class \code{MultiAssayExperiment}, the factor of sample classes must be stored
  in the DataFrame accessible by the \code{colData} function with column name \code{"class"}.
}
\value{
  Either a factor vector of predicted classes, a matrix of scores for each class, or a table of
  both the class labels and class scores, depending on the setting of \code{returnType}.
}
\author{Dario Strbenac}
\examples{
  if(require(PoiClaClu))
  {
    readCounts <- CountDataSet(n = 100, p = 1000, 2, 5, 0.1)
    # Rows are for features, columns are for samples.
    trainData <- t(readCounts[['x']])
    classes <- factor(paste("Class", readCounts[['y']]))
    testData <- t(readCounts[['xte']])
    storage.mode(trainData) <- storage.mode(testData) <- "integer"
    classified <- classifyInterface(trainData, classes, testData)
    
    setNames(table(paste("Class", readCounts[["yte"]]) == classified), c("Incorrect", "Correct"))
  }
}
