#' Ranking of Differentially Abundant Features
#' 
#' Uses an ordinary t-test if the data set has two classes or one-way ANOVA if
#' the data set has three or more classes to select differentially expressed
#' features.
#' 
#' This ranking method looks for changes in means and uses
#' \code{\link[genefilter]{rowttests}} to rank the features if there are two
#' classes or \code{\link[genefilter]{rowFtests}} if there are three or more
#' classes.
#' 
#' @aliases differentMeansRanking differentMeansRanking,matrix-method
#' differentMeansRanking,DataFrame-method
#' differentMeansRanking,MultiAssayExperiment-method
#' @param measurements Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix}, the rows are features, and the columns are samples.
#' @param classes A vector of class labels of class \code{\link{factor}} of the
#' same length as the number of samples in \code{measurements} if it is a
#' \code{\link{matrix}} (i.e. number of columns) or a \code{\link{DataFrame}}
#' (i.e. number of rows) or a character vector of length 1 containing the
#' column name in \code{measurements} if it is a \code{\link{DataFrame}} or the
#' column name in \code{colData(measurements)} if \code{measurements} is a
#' \code{\link{MultiAssayExperiment}}. If a column name, that column will be
#' removed before training.
#' @param targets Names of data tables to be combined into a single table and
#' used in the analysis.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return A vector or data frame (if \code{MultiAssayExperiment} input) of
#' features, from the most promising features in the first position to the
#' least promising feature in the last position.
#' @author Dario Strbenac
#' @examples
#' 
#'     # Genes 76 to 100 have differential expression.
#'     genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
#'     genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
#'                                  c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix))
#'     rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
#'     
#'     ranked <- differentMeansRanking(genesMatrix, classes)
#'     head(ranked)
#' 
#' @export
setGeneric("differentMeansRanking", function(measurements, ...)
           standardGeneric("differentMeansRanking"))

# Matrix of numeric measurements.
setMethod("differentMeansRanking", "matrix", function(measurements, classes, ...)
{
  differentMeansRanking(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

# DataFrame of numeric measurements, likely created by runTests or runTest.
setMethod("differentMeansRanking", "DataFrame",
          function(measurements, classes, verbose = 3)
{
  if(!requireNamespace("genefilter", quietly = TRUE))
    stop("The package 'genefilter' could not be found. Please install it.")

  measurementsMatrix <- t(as.matrix(measurements))
  
  if(length(levels(classes)) == 2)
    pValues <- genefilter::rowttests(measurementsMatrix, classes)[, "p.value"]
  else
    pValues <- genefilter::rowFtests(measurementsMatrix, classes)[, "p.value"]
  
  if(!is.null(S4Vectors::mcols(measurements)))
    S4Vectors::mcols(measurements)[order(pValues), ]
  else
    colnames(measurements)[order(pValues)]
})

# One or more omics data sets, possibly with clinical data.
setMethod("differentMeansRanking", "MultiAssayExperiment", 
          function(measurements, targets = NULL, classes, ...)
{
  if(is.null(targets))
    stop("'targets' must be specified but was not.")
  if(length(setdiff(targets, names(measurements))))
    stop("Some values of 'targets' are not names of 'measurements' but all must be.")                            
            
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  differentMeansRanking(measurements, classes, ...)
})