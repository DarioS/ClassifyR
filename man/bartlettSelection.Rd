\name{bartlettSelection}
\alias{bartlettSelection}
\alias{bartlettSelection,matrix-method}
\alias{bartlettSelection,DataFrame-method}
\alias{bartlettSelection,MultiAssayExperiment-method}
\title{Selection of Differential Variability with Bartlett Statistic}
\description{Ranks features by largest Bartlett statistic and chooses
  the features which have best resubstitution performance.
}
\usage{
  \S4method{bartlettSelection}{matrix}(measurements, classes, ...)
  \S4method{bartlettSelection}{DataFrame}(measurements, classes, ...)
  \S4method{bartlettSelection}{MultiAssayExperiment}(measurements, targets, ...)
  (measurements, targets = names(measurements), datasetName,
  trainParams, predictParams, resubstituteParams,
  selectionName = "Bartlett Test", verbose = 3)
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
    \code{\link{MultiAssayExperiment}} containing the training data. For a \code{matrix},
    the rows are features, and the columns are samples.}
  \item{targets}{The names of the data tables to be used. \code{"clinical"} is also a valid value and specifies that numeric variables from the clinical data table will be used.}    
  \item{classes}{Either a vector of class labels of class \code{\link{factor}} or if
                 the measurements are of class \code{DataFrame} a character vector of
                 length 1 containing the column name in \code{measurement} is also permitted.
                 Not used if \code{measurements} is a \code{MultiAssayExperiment} object.}
  \item{targets}{Specifies the data tables to consider If \code{measurements} is a
                \code{MultiAssayExperiment}. A vector of data table names that is a subset
                of \code{names(measurements)}.}
  \item{...}{Unused variables by the three top-level methods passed to the internal method
             which generates the plot(s).}
  \item{datasetName}{A name for the dataset used. Stored in the result.}
  \item{trainParams}{A container of class \code{\link{TrainParams}} describing the
    classifier to use for training.}
  \item{predictParams}{A container of class \code{\link{PredictParams}} describing how
    prediction is to be done.}
  \item{resubstituteParams}{An object of class \code{\link{ResubstituteParams}}
    describing the performance measure to consider and the numbers of
    top features to try for resubstitution classification.}                               
  \item{selectionName}{A name to identify this selection method by. Stored in the result.}
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
    This function only prints progress messages if the value is 3.}
}
\details{
  The calculation of the test statistic is performed by the \code{\link{bartlett.test}}
  function from the \code{\link{stats}} package.
  
  Data tables which consist entirely of non-numeric data must not be specified by the
  \code{targets} parameter.
}
\value{
  An object of class \code{\link{SelectResult}} or a list of such objects, if the classifier
  which was used for determining resubstitution error rate made a number of prediction varieties.
}
\author{Dario Strbenac}

\examples{
  # Samples in one class with differential variability to other class.
  # First 20 genes are DV.
  genesRNAmatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 1)))
  moreVariable <- sapply(1:25, function(sample) rnorm(20, 9, 5))
  genesRNAmatrix <- cbind(genesRNAmatrix, rbind(moreVariable,
                          sapply(1:25, function(sample) rnorm(80, 9, 1))))
  colnames(genesRNAmatrix) <- paste("Sample", 1:50)
  rownames(genesRNAmatrix) <- paste("Gene", 1:100)
  genesSNPmatrix <- matrix(sample(c("None", "Missense"), 250, replace = TRUE),
                           ncol = 50)
  colnames(genesSNPmatrix) <- paste("Sample", 1:50)
  rownames(genesSNPmatrix) <- paste("Gene", 1:5)
  classes <- factor(rep(c("Poor", "Good"), each = 25))
  names(classes) <- paste("Sample", 1:50)
  genesDataset <- MultiAssayExperiment(list(RNA = genesRNAmatrix, SNP = genesSNPmatrix),
                                       colData = DataFrame(class = classes))
    
  trainIDs <- paste("Sample", c(1:20, 26:45))
  genesDataset <- subtractFromLocation(genesDataset, training = trainIDs,
                                        targets = "RNA") # Exclude SNP data.
                                         
  resubParams <- ResubstituteParams(nFeatures = seq(10, 100, 10),
                                    performanceType = "balanced", better = "lower")
  bartlettSelection(genesDataset, datasetName = "Example", targets = "RNA",
                    trainParams = TrainParams(fisherDiscriminant, FALSE, TRUE),
                    predictParams = PredictParams(function(){}, FALSE,
                    getClasses = function(result) result),
                    resubstituteParams = resubParams)
}
