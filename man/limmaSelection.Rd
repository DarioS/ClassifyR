\name{limmaSelection}
\alias{limmaSelection}
\alias{limmaSelection,matrix-method}
\alias{limmaSelection,DataFrame-method}
\alias{limmaSelection,MultiAssayExperiment-method}
\title{Selection of Differentially Abundant Features}
\description{Uses a moderated t-test with empirical Bayes shrinkage to select
             differentially expressed features.
       }
\usage{
  \S4method{limmaSelection}{matrix}(measurements, classes, ...)
  \S4method{limmaSelection}{DataFrame}(measurements, classes, datasetName,
               trainParams, predictParams, resubstituteParams, ...,
               selectionName = "Moderated t-test", verbose = 3)
  \S4method{limmaSelection}{MultiAssayExperiment}(measurements, targets = NULL, ...)
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
                      \code{\link{MultiAssayExperiment}} containing the training data.
                      For a \code{matrix}, the rows are features, and the columns are samples.}
  \item{classes}{A vector of class labels of class \code{\link{factor}} of the same length
                 as the number of samples in \code{measurements}. Not used if \code{measurements}
                 is a \code{MultiAssayExperiment} object.}
  \item{targets}{Names of data tables to be combined into a single table and
                 used in the analysis.}                  
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method or optional settings that are
            passed to \code{\link[limma]{lmFit}}.}
  \item{datasetName}{A name for the data set used. Stored in the result.}                    
  \item{trainParams}{A container of class \code{\link{TrainParams}} describing the
                     classifier to use for training.}
  \item{predictParams}{A container of class \code{\link{PredictParams}} describing how
                       prediction is to be done.}
  \item{resubstituteParams}{An object of class \code{\link{ResubstituteParams}}
                            describing the performance measure to consider and the numbers of
                            top features to try for resubstitution classification.}                     
  \item{selectionName}{A name to identify this selection method by. Stored in the result.}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  This selection method looks for changes in means and uses a moderated t-test.
}
\value{
  An object of class \code{\link{SelectResult}} or a list of such objects, if the classifier
  which was used for determining the specified performance metric made a number of prediction varieties.
}

\author{Dario Strbenac}
\references{
Limma: linear models for microarray data, Gordon Smyth, 2005, In: Bioinformatics and Computational Biology Solutions using R and Bioconductor, Springer, New York, pages 397-420.
}

\examples{
  #if(require(sparsediscrim))
  #{
    # Genes 76 to 100 have differential expression.
    genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
                                 c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix))
    rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
    
    resubstituteParams <- ResubstituteParams(nFeatures = seq(10, 100, 10),
                            performanceType = "balanced error", better = "lower")
    selected <- limmaSelection(genesMatrix, classes, "Example",
                               trainParams = TrainParams(), predictParams = PredictParams(),
                               resubstituteParams = resubstituteParams)
                               
    selected@chosenFeatures
  #}
}
