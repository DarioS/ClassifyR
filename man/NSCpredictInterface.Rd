\name{NSCpredictInterface}
\alias{NSCpredictInterface}
\alias{NSCpredictInterface,pamrtrained,matrix-method}
\alias{NSCpredictInterface,pamrtrained,DataFrame-method}
\alias{NSCpredictInterface,pamrtrained,MultiAssayExperiment-method}

\title{Interface for \code{pamr.predict} Function from \code{pamr} CRAN Package}
\description{Restructures variables from ClassifyR framework to be compatible with \code{\link[pamr]{pamr.predict}} definition.}
\usage{
  \S4method{NSCpredictInterface}{pamrtrained,matrix}(trained, test, ...)
  \S4method{NSCpredictInterface}{pamrtrained,DataFrame}(trained, test, classes = NULL, ...,
                                 returnType = c("both", "class", "score"),
                                 classifierName = "Nearest Shrunken Centroids", verbose = 3)  
  \S4method{NSCpredictInterface}{pamrtrained,MultiAssayExperiment}(trained, test, targets = names(test), ...)
}
\arguments{
  \item{trained}{An object of class \code{pamrtrained}.}             
  \item{test}{An object of the same class as \code{measurements} with no samples in common with
              \code{measurements} used in the training stage and the same number of features as it.
              Also, if a \code{DataFrame}, the \code{class} column must be absent.}
  \item{classes}{Either NULL or a character vector of length 1, specifying the column name to remove.}            
  \item{targets}{If \code{test} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 numeric variables from the clinical data table will be used.}   
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method or optional settings that are
            passed to \code{\link[pamr]{pamr.predict}}.}
  \item{returnType}{Default: \code{"both"}. Either \code{"class"}, \code{"score"} or \code{"both"}.
                    Sets the return value from the prediction to either a vector of class labels,
                    score for a sample belonging to the second class, as determined by the
                    factor levels, or both labels and scores in a \code{data.frame}.}
  \item{classifierName}{Default: Nearest Shrunken Centroids. Useful for automated plot annotation by
                        plotting functions within this package.}                    
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  This function is an interface between the ClassifyR framework and \code{\link[pamr]{pamr.predict}}.
  It selects the highest threshold that gives the minimum error rate in the training data.
}
\value{
  Either a factor vector of predicted classes, a matrix of scores for each class, or a table of
  both the class labels and class scores, depending on the setting of \code{returnType}.
}
\author{Dario Strbenac}
\seealso{\code{\link[pamr]{pamr.predict}} for the function that was interfaced to.}
\examples{
  if(require(pamr))
  {
    # Samples in one class with differential expression to other class.
    genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(100, 9, 1)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn)
                                 c(rnorm(75, 9, 1), rnorm(25, 14, 1))))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    
    fit <- NSCtrainInterface(genesMatrix[, c(1:20, 26:45)], classes[c(1:20, 26:45)])
    NSCpredictInterface(fit, genesMatrix[, c(21:25, 46:50)])
  }
}
