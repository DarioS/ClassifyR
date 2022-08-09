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
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data.  For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
#' If of type \code{\link{DataFrame}} or \code{\link{MultiAssayExperiment}}, the data set is subset
#' to only those features of type \code{numeric}.
#' @param classesTrain A vector of class labels of class \code{\link{factor}} of the
#' same length as the number of samples in \code{measurementsTrain} if it is a
#' \code{\link{matrix}} or a \code{\link{DataFrame}} or a character vector of length 1
#' containing the column name in \code{measurementsTrain} if it is a \code{\link{DataFrame}} or the
#' column name in \code{colData(measurementsTrain)} if \code{measurementsTrain} is a
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
#' @return A vector of feature indices, from the most promising features in the
#' first position to the least promising feature in the last position.
#' @author Dario Strbenac
#' @examples
#' 
#'     # Genes 76 to 100 have differential expression.
#'     genesMatrix <- sapply(1:100, function(sample) rnorm(25, 9, 0.3))
#'     genesMatrix <- rbind(genesMatrix, t(sapply(1:25, function(sample)
#'                                       c(rnorm(75, 9, 0.3), rnorm(25, 14, 0.3)))))
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     rownames(genesMatrix) <- paste("Sample", 1:50)
#'     colnames(genesMatrix) <- paste("Gene", 1:100)
#'     
#'     ranked <- differentMeansRanking(genesMatrix, classes)
#'     head(ranked)
#' @usage NULL
#' @export
#' @importFrom genefilter rowttests rowFtests
setGeneric("differentMeansRanking", function(measurementsTrain, ...)
           standardGeneric("differentMeansRanking"))

# Matrix of numeric measurements.
#' @rdname differentMeansRanking
#' @export
setMethod("differentMeansRanking", "matrix", function(measurementsTrain, classesTrain, ...)
{
  differentMeansRanking(S4Vectors::DataFrame(measurementsTrain, check.names = FALSE), classesTrain, ...)
})

# DataFrame of numeric measurements, likely created by runTests or runTest.
#' @rdname differentMeansRanking
#' @export
setMethod("differentMeansRanking", "DataFrame",
          function(measurementsTrain, classesTrain, verbose = 3)
{
  if(!requireNamespace("genefilter", quietly = TRUE))
    stop("The package 'genefilter' could not be found. Please install it.")

  splitDataset <- .splitDataAndOutcome(measurementsTrain, classesTrain)
  classesTrain <- splitDataset[["outcome"]]
  # Data is required to be in traditional bioinformatics format - features in rows
  # and samples in columns and also must be a matrix, not another kind of rectangular data.  
  
  pValues <- NULL
  
  categ <- sapply(splitDataset[["measurements"]], class) %in% c("character", "factor")
  if(any(categ)){
     pValues[categ] <- sapply(which(categ), function(x){
      chisq.test(splitDataset[["measurements"]][,x], classesTrain)$p.value
    })
  }
  

  
  
  if(any(!categ)){
    measurementsMatrix <- t(as.matrix(splitDataset[["measurements"]][,!categ, drop = FALSE]))
  if(length(levels(classesTrain)) == 2)
  {
    if(verbose == 3)
      message("Ranking features based on t-statistic.")
    pValues[!categ] <- genefilter::rowttests(measurementsMatrix, classesTrain)[, "p.value"]
  } else {
    if(verbose == 3)
      message("Ranking features based on F-statistic.")
    pValues[!categ]  <- genefilter::rowFtests(measurementsMatrix, classesTrain)[, "p.value"]
  }
  }
  
  order(pValues) # From smallest to largest.
})

# One or more omics data sets, possibly with sample information data.
#' @rdname differentMeansRanking
#' @export
setMethod("differentMeansRanking", "MultiAssayExperiment", 
          function(measurementsTrain, targets = NULL, classesTrain, ...)
{
  if(is.null(targets))
    stop("'targets' must be specified but was not.")
  if(length(setdiff(targets, names(measurementsTrain))))
    stop("Some values of 'targets' are not names of 'measurements' but all must be.")                            
            
  tablesAndClasses <- .MAEtoWideTable(measurementsTrain, targets, classesTrain)
  measurementsTrain <- tablesAndClasses[["dataTable"]]
  classesTrain <- tablesAndClasses[["outcome"]]
  differentMeansRanking(measurementsTrain, classesTrain, ...)
})