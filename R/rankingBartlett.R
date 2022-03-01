#' Ranking of Differential Variability with Bartlett Statistic
#' 
#' Ranks all features from largest Bartlett statistic to smallest.
#' 
#' The calculation of the test statistic is performed by the
#' \code{\link{bartlett.test}} function from the \code{\link{stats}} package.
#' 
#' Data tables which consist entirely of non-numeric data cannot be ranked.
#' 
#' @aliases bartlettRanking bartlettRanking,matrix-method
#' bartlettRanking,DataFrame-method bartlettRanking,MultiAssayExperiment-method
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
#' @param targets If \code{measurementsTrain} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"sampleInfo"} is also a valid value
#' and specifies that numeric variables from the sample information table will be
#' used.
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
#'   # Samples in one class with differential variability to other class.
#'   # First 20 genes are DV.
#'   genesRNAmatrix <- sapply(1:20, function(sample) c(rnorm(25, 9, 1), rnorm(25, 9, 5)))
#'   genesRNAmatrix <- cbind(genesRNAmatrix, sapply(1:80, function(sample) rnorm(50, 9, 1)))
#'   rownames(genesRNAmatrix) <- paste("Sample", 1:50)
#'   colnames(genesRNAmatrix) <- paste("Gene", 1:100)
#'   genesSNPmatrix <- matrix(sample(c("None", "Missense"), 250, replace = TRUE), nrow = 50)
#'   rownames(genesSNPmatrix) <- paste("Sample", 1:50)
#'   colnames(genesSNPmatrix) <- paste("Gene", 1:5)
#'   classes <- factor(rep(c("Poor", "Good"), each = 25))
#'   names(classes) <- paste("Sample", 1:50)
#'   
#'   
#'   genesDataset <- MultiAssayExperiment(list(RNA = t(genesRNAmatrix), SNP = t(genesSNPmatrix)),
#'                                        colData = DataFrame(class = classes))
#'                                          
#'   bartlettRanking(genesDataset, targets = "RNA", classesTrain = "class")
#'
#' @import stats methods
#' @export
setGeneric("bartlettRanking", function(measurementsTrain, ...)
standardGeneric("bartlettRanking"))

setMethod("bartlettRanking", "matrix", # Matrix of numeric measurements.
function(measurementsTrain, classesTrain, ...)
{
  bartlettRanking(DataFrame(measurementsTrain, check.names = FALSE), classesTrain, ...)
})

setMethod("bartlettRanking", "DataFrame", # Sample information data or one of the other inputs, transformed.
          function(measurementsTrain, classesTrain, verbose = 3)
{
  splitDataset <- .splitDataAndOutcomes(measurementsTrain, classesTrain)
  measurementsTrain <- splitDataset[["measurements"]]
  classesTrain <- splitDataset[["outcomes"]]
  
  if(verbose == 3)
    message("Ranking features based on Bartlett statistic.")
  
  pValues <- apply(measurementsTrain, 2, function(featureColumn)
    stats::bartlett.test(featureColumn, classesTrain)[["p.value"]])
  
  if(!is.null(S4Vectors::mcols(measurementsTrain)))
    S4Vectors::mcols(measurementsTrain)[order(pValues), ]
  else
    colnames(measurementsTrain)[order(pValues)]
})

# One or more omics data sets, possibly with sample information data.
setMethod("bartlettRanking", "MultiAssayExperiment",
          function(measurementsTrain, targets = names(measurementsTrain), classesTrain, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurementsTrain, targets, classesTrain)
  measurementsTrain <- tablesAndClasses[["dataTable"]]
  classesTrain <- tablesAndClasses[["outcomes"]]
  
  if(ncol(measurementsTrain) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    bartlettRanking(measurementsTrain, classesTrain, ...)
})