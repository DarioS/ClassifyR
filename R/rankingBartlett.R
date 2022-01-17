#' Ranking of Differential Variability with Bartlett Statistic
#' 
#' Ranks allfeatures from largest Bartlett statistic to smallest.
#' 
#' The calculation of the test statistic is performed by the
#' \code{\link{bartlett.test}} function from the \code{\link{stats}} package.
#' 
#' Data tables which consist entirely of non-numeric data cannot be ranked.
#' 
#' @aliases bartlettRanking bartlettRanking,matrix-method
#' bartlettRanking,DataFrame-method bartlettRanking,MultiAssayExperiment-method
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
#' removed before ranking.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"clinical"} is also a valid value
#' and specifies that numeric variables from the clinical data table will be
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
#'   genesRNAmatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 1)))
#'   moreVariable <- sapply(1:25, function(sample) rnorm(20, 9, 5))
#'   genesRNAmatrix <- cbind(genesRNAmatrix, rbind(moreVariable,
#'                           sapply(1:25, function(sample) rnorm(80, 9, 1))))
#'   colnames(genesRNAmatrix) <- paste("Sample", 1:50)
#'   rownames(genesRNAmatrix) <- paste("Gene", 1:100)
#'   genesSNPmatrix <- matrix(sample(c("None", "Missense"), 250, replace = TRUE),
#'                            ncol = 50)
#'   colnames(genesSNPmatrix) <- paste("Sample", 1:50)
#'   rownames(genesSNPmatrix) <- paste("Gene", 1:5)
#'   classes <- factor(rep(c("Poor", "Good"), each = 25))
#'   names(classes) <- paste("Sample", 1:50)
#'   genesDataset <- MultiAssayExperiment(list(RNA = genesRNAmatrix, SNP = genesSNPmatrix),
#'                                        colData = DataFrame(class = classes))
#'   # Wait for update to MultiAssayExperiment wideFormat function.  
#'   trainIDs <- paste("Sample", c(1:20, 26:45))
#'   genesDataset <- subtractFromLocation(genesDataset, training = trainIDs,
#'                                        targets = "RNA") # Exclude SNP data.
#'                                          
#'   bartlettRanking(genesDataset, targets = "RNA")
#'
#' @import stats methods
#' @export
setGeneric("bartlettRanking", function(measurements, ...)
standardGeneric("bartlettRanking"))

setMethod("bartlettRanking", "matrix", # Matrix of numeric measurements.
function(measurements, classes, ...)
{
  bartlettRanking(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("bartlettRanking", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  classes <- splitDataset[["classes"]]
  isNumeric <- sapply(measurements, is.numeric)
  measurements <- measurements[, isNumeric, drop = FALSE]
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  
  if(verbose == 3)
    message("Ranking features based on Bartlett statistic.")
  
  pValues <- apply(measurements, 2, function(featureColumn)
    stats::bartlett.test(featureColumn, classes)[["p.value"]])
  
  if(!is.null(S4Vectors::mcols(measurements)))
    S4Vectors::mcols(measurements)[order(pValues), ]
  else
    colnames(measurements)[order(pValues)]
})

# One or more omics data sets, possibly with clinical data.
setMethod("bartlettRanking", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), classes, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  
  if(ncol(measurements) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    bartlettRanking(measurements, classes, ...)
})