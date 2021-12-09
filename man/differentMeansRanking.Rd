\name{differentMeansRanking}
\alias{differentMeansRanking}
\alias{differentMeansRanking,matrix-method}
\alias{differentMeansRanking,DataFrame-method}
\alias{differentMeansRanking,MultiAssayExperiment-method}
\title{Ranking of Differentially Abundant Features}
\description{Uses an ordinary t-test if the data set has two classes or one-way ANOVA if the data set has three or more classes to select differentially expressed features.
       }
\usage{
  \S4method{differentMeansRanking}{matrix}(measurements, classes, ...)
  \S4method{differentMeansRanking}{DataFrame}(measurements, classes, verbose = 3)
  \S4method{differentMeansRanking}{MultiAssayExperiment}(measurements, targets = NULL, ...)
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
             are passed into and used by the \code{DataFrame} method.}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  This ranking method looks for changes in means and uses \code{\link[genefilter]{rowttests}} to rank the features if there are two classes or \code{\link[genefilter]{rowFtests}} if there are three or more classes.
}
\value{
  A vector or data frame (if \code{MultiAssayExperiment} input) of features, from the most promising features in the first position to the least promising feature in the last position.
}
\author{Dario Strbenac}
\examples{
    # Genes 76 to 100 have differential expression.
    genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
                                 c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix))
    rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
    
    ranked <- differentMeansRanking(genesMatrix, classes)
    head(ranked)
}
