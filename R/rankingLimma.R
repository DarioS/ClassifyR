#' Ranking of Differentially Abundant Features
#' 
#' Uses a moderated F-test with empirical Bayes shrinkage to rank
#' differentially expressed features based on differences of means. This means
#' it works when there are three or more classes.
#' 
#' This ranking method looks for changes in means and uses a moderated F-test
#' to do so.
#' 
#' @aliases limmaRanking limmaRanking,matrix-method
#' limmaRanking,DataFrame-method limmaRanking,MultiAssayExperiment-method
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
#' If of type \code{\link{DataFrame}} or \code{\link{MultiAssayExperiment}}, the data set is subset
#' to only those features of type \code{numeric}.
#' @param classesTrain A vector of class labels of class \code{\link{factor}} of the
#' same length as the number of samples in \code{measurements}.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method or optional settings that are passed to
#' \code{\link[limma]{lmFit}}.
#' @param targets Names of data tables to be combined into a single table and
#' used in the analysis.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return A vector or data frame (if \code{MultiAssayExperiment} input) of
#' features, from the most promising features in the first position to the
#' least promising feature in the last position.
#' @author Dario Strbenac
#' @references Limma: linear models for microarray data, Gordon Smyth, 2005,
#' In: Bioinformatics and Computational Biology Solutions using R and
#' Bioconductor, Springer, New York, pages 397-420.
#' @examples
#' 
#'   #if(require(sparsediscrim))
#'   #{
#'     # Genes 76 to 100 have differential expression.
#'     genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
#'     genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
#'                                  c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix))
#'     rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
#'     
#'     ranked <- limmaRanking(genesMatrix, classes)
#'     head(ranked)                           
#'   #}
#' 
#' @export
setGeneric("limmaRanking", function(measurements, ...)
           standardGeneric("limmaRanking"))

# Matrix of numeric measurements.
setMethod("limmaRanking", "matrix", function(measurements, classesTrain, ...)
{
  limmaRanking(DataFrame(measurementsTrain, check.names = FALSE), classesTrain, ...)
})

# DataFrame of numeric measurements, likely created by runTests or runTest.
setMethod("limmaRanking", "DataFrame",
          function(measurementsTrain, classesTrain, ..., verbose = 3)
{
  if(!requireNamespace("limma", quietly = TRUE))
    stop("The package 'limma' could not be found. Please install it.")

  fitParams <- list(t(as.matrix(measurementsTrain)), model.matrix(~ classesTrain))
  if(!missing(...))
    fitParams <- append(fitParams, ...)
  linearModel <- do.call(limma::lmFit, fitParams)
  linearModel <- limma::eBayes(linearModel)
  linearModel <- linearModel[, -1] # Get rid of intercept.
  
  if(!is.null(S4Vectors::mcols(measurementsTrain)))
    S4Vectors::mcols(measurementsTrain)[order(linearModel[["F.p.value"]]), ]
  else
    colnames(measurementsTrain)[order(linearModel[["F.p.value"]])]
})

# One or more omics data sets, possibly with sample information data.
setMethod("limmaRanking", "MultiAssayExperiment", 
          function(measurementsTrain, targets = NULL, classesTrain, ...)
{
  if(is.null(targets))
    stop("'targets' must be specified but was not.")
  if(length(setdiff(targets, names(measurementsTrain))))
    stop("Some values of 'targets' are not names of 'measurementsTrain' but all must be.")                            
            
  tablesAndClasses <- .MAEtoWideTable(measurementsTrain, targets, classesTrain)
  measurementsTrain <- tablesAndClasses[["dataTable"]]
  classesTrain <- tablesAndClasses[["outcomes"]]
  limmaRanking(measurementsTrain, classesTrain, ...)
})