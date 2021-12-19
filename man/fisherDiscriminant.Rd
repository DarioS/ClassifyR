\name{fisherDiscriminant}
\alias{fisherDiscriminant}
\alias{fisherDiscriminant,matrix-method}
\alias{fisherDiscriminant,DataFrame-method}
\alias{fisherDiscriminant,MultiAssayExperiment-method}
\title{Classification Using Fisher's LDA}
\description{Finds the decision boundary using the training set, and gives predictions
for the test set.
}
\usage{
  \S4method{fisherDiscriminant}{matrix}(measurements, classes, test, ...)
  \S4method{fisherDiscriminant}{DataFrame}(measurements, classes, test, returnType = c("both", "class", "score"),
                   verbose = 3)
  \S4method{fisherDiscriminant}{MultiAssayExperiment}(measurements, test, targets = names(measurements), classes, ...)
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
                       \code{\link{MultiAssayExperiment}} containing the training data.
                       For a \code{matrix}, the rows are features, and the columns are samples.}
  \item{classes}{A vector of class labels of class \code{\link{factor}} of the same length
                 as the number of samples in \code{measurements} if it is a \code{\link{matrix}} (i.e. number of columns) or a \code{\link{DataFrame}} (i.e. number of rows) or a character vector
                 of length 1 containing the column name in \code{measurements} if it is a \code{\link{DataFrame}} or the column name in \code{colData(measurements)} if \code{measurements} is
a \code{\link{MultiAssayExperiment}}. If a column name, that column will be removed before training.}
  \item{test}{An object of the same class as \code{measurements} with no samples in common with
              \code{measurements} and the same number of features as it.}
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 integer variables from the clinical data table will be used.}              
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method.}
  \item{returnType}{Default: \code{"both"}. Either \code{"class"}, \code{"score"}, or \code{"both"}.
                    Sets the return value from the prediction to either a vector of class labels,
                    score for a sample belonging to the second class, as determined by the factor levels,
                    or both labels and scores in a \code{data.frame}.}                    
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  Unlike ordinary LDA, Fisher's version does not have assumptions about the normality
  of the features.
  
  Data tables which consist entirely of non-numeric data cannot be analysed. If \code{measurements}
  is an object of class \code{MultiAssayExperiment}, the factor of sample classes must be stored
  in the DataFrame accessible by the \code{colData} function with column name \code{"class"}.
}
\value{
  A vector or \code{data.frame} of class prediction information, as long as the number of samples in the test data.
}
\author{Dario Strbenac}

\examples{
  trainMatrix <- matrix(rnorm(1000, 8, 2), ncol = 10)
  classes <- factor(rep(c("Poor", "Good"), each = 5))
  
  # Make first 30 genes increased in value for poor samples.
  trainMatrix[1:30, 1:5] <- trainMatrix[1:30, 1:5] + 5
  
  testMatrix <- matrix(rnorm(1000, 8, 2), ncol = 10)
  
  # Make first 30 genes increased in value for sixth to tenth samples.
  testMatrix[1:30, 6:10] <- testMatrix[1:30, 6:10] + 5
  
  fisherDiscriminant(trainMatrix, classes, testMatrix)
}
