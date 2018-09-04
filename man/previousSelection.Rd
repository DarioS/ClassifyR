\name{previousSelection}
\alias{previousSelection}
\alias{previousSelection,matrix-method}
\alias{previousSelection,DataFrame-method}
\alias{previousSelection,MultiAssayExperiment-method}
\title{Automated Selection of Previously Selected Features}
\description{Uses the feature selection of the same cross-validation iteration of a previous classification
for the current classification task.
}
\usage{
  \S4method{previousSelection}{matrix}(measurements, ...)
  \S4method{previousSelection}{DataFrame}(measurements, classes, datasetName,
                  classifyResult, minimumOverlapPercent = 80,
                        selectionName = "Previous Selection", .iteration, verbose = 3)
  \S4method{previousSelection}{MultiAssayExperiment}(measurements, ...)
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
                      \code{\link{MultiAssayExperiment}} containing the training data.
                      For a \code{matrix}, the rows are features, and the columns are samples.}
  \item{classes}{Do not specify this variable. It is ignored and only used to create consistency of formal
                 parameters with other feature selection methods.}                      
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method.}  
  \item{datasetName}{A name for the data set used. Stored in the result.}
  \item{classifyResult}{An existing classification result from which to take the feature selections from.}
  \item{minimumOverlapPercent}{If at least this many selected features can't be identified in the
                               current data set, then the selection stops with an error.}
  \item{selectionName}{A name to identify this selection method by. Stored in the result.}
  \item{.iteration}{Do not specify this variable. It is set by \code{\link{runTests}} if this function
                    is being repeatedly called by \code{runTests}.}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
    This function only prints progress messages if the value is 3.}
}
\value{
  An object of class \code{\link{SelectResult}}.
}
\author{Dario Strbenac}

\examples{
  #if(require(sparsediscrim))
  #{
    # Genes 76 to 100 have differential expression.
    genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
                                 c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
    colnames(genesMatrix) <- paste("Sample", 1:50)
    rownames(genesMatrix) <- paste("Gene", 1:100)                                 
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    resubstitute <- ResubstituteParams(nFeatures = seq(10, 100, 10), 
                        performanceType = "error", better = "lower")
    result <- runTests(genesMatrix, classes, datasetName = "Example",
                       classificationName = "Differential Expression",
                       permutations = 2, fold = 2,
                       params = list(SelectParams(), TrainParams(), PredictParams()))
                       
    # Genes 50 to 74 have differential expression in new data set.
    newDataset <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
    newDataset <- cbind(newDataset, rbind(sapply(1:25, function(sample) rnorm(49, 9, 2)),
                                          sapply(1:25, function(sample) rnorm(25, 14, 2)),
                                          sapply(1:25, function(sample) rnorm(26, 9, 2))))
                                          
    rownames(newDataset) <- rownames(genesMatrix)
    colnames(newDataset) <- colnames(genesMatrix)
                                           
    newerResult <- runTests(newDataset, classes, datasetName = "Latest Data Set",
                            classificationName = "Differential Expression",
                            permutations = 2, fold = 2,
                       params = list(SelectParams(previousSelection,
                                     intermediate = ".iteration",
                                     classifyResult = result),
                       TrainParams(DLDAtrainInterface),
                       PredictParams(DLDApredictInterface,
                                     getClasses = function(result) result[["class"]])))
                                     
    # However, only genes 76 to 100 are chosen, because the feature selections are
    # carried over from the first cross-validated classification.
    features(newerResult)
  #}  
}
