#' Ranking of Differential Variability with coxph Statistic
#' 
#' Ranks allfeatures from largest coxph statistic to smallest.
#' 
#' The calculation of the test statistic is performed by the
#' \code{\link{coxph.test}} function from the \code{\link{stats}} package.
#' 
#' Data tables which consist entirely of non-numeric data cannot be ranked.
#' 
#' @aliases coxphRanking coxphRanking,matrix-method
#' coxphRanking,DataFrame-method coxphRanking,MultiAssayExperiment-method
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
#' @import survival coxph
#' @export
setGeneric("coxphRanking", function(measurements, ...)
standardGeneric("coxphRanking"))

setMethod("coxphRanking", "matrix", # Matrix of numeric measurements.
function(measurements, classes, ...)
{
  coxphRanking(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("coxphRanking", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, verbose = 3)
{
  splitDataset <- .splitDataAndClasses(measurements, classes)
  measurements <- splitDataset[["measurements"]]
  classes <- splitDataset[["classes"]]

  pValues <- apply(measurements, 2, function(featureColumn){
    fit <- survival::coxph(classes ~ featureColumn)
    s <- summary(fit)
    s$waldtest["pvalue"]
  })
  
  if(!is.null(S4Vectors::mcols(measurements)))
    S4Vectors::mcols(measurements)[order(pValues), ]
  else
    colnames(measurements)[order(pValues)]
})

# One or more omics data sets, possibly with clinical data.
setMethod("coxphRanking", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), classes, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes)
  measurements <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]
  
  if(ncol(measurements) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    coxphRanking(measurements, classes, ...)
})