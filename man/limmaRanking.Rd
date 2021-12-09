\name{limmaRanking}
\alias{limmaRanking}
\alias{limmaRanking,matrix-method}
\alias{limmaRanking,DataFrame-method}
\alias{limmaRanking,MultiAssayExperiment-method}
\title{Ranking of Differentially Abundant Features}
\description{Uses a moderated F-test with empirical Bayes shrinkage to rank
             differentially expressed features based on differences of means.
             This means it works when there are three or more classes.
       }
\usage{
  \S4method{limmaRanking}{matrix}(measurements, classes, ...)
  \S4method{limmaRanking}{DataFrame}(measurements, classes, ..., verbose = 3)
  \S4method{limmaRanking}{MultiAssayExperiment}(measurements, targets = NULL, ...)
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
                      \code{\link{MultiAssayExperiment}} containing the training data.
                      For a \code{matrix}, the rows are features, and the columns are samples.}
  \item{classes}{A vector of class labels of class \code{\link{factor}} of the same length
                 as the number of samples in \code{measurements}. Not used if \code{measurements}
                 is a \code{MultiAssayExperiment} object.}
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method or optional settings that are
            passed to \code{\link[limma]{lmFit}}.}
  \item{targets}{Names of data tables to be combined into a single table and
                 used in the analysis.}            
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  This ranking method looks for changes in means and uses a moderated F-test to do so.
}
\value{
  A vector or data frame (if \code{MultiAssayExperiment} input) of features, from the most promising features in the first position to the least promising feature in the last position.
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
    
    ranked <- limmaRanking(genesMatrix, classes)
    head(ranked)                           
  #}
}
