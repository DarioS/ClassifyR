#' Ranking of Differential Variability with coxph Statistic
#' 
#' Ranks all features from largest coxph statistic to smallest.
#' 
#' The calculation of the test statistic is performed by the
#' \code{\link{coxph.test}} function from the \code{\link{stats}} package.
#' 
#' Data tables which consist entirely of non-numeric data cannot be ranked.
#' 
#' @aliases coxphRanking coxphRanking,matrix-method
#' coxphRanking,DataFrame-method coxphRanking,MultiAssayExperiment-method
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
#' @param survivalTrain A tabular data type of survival information of the
#' same number of rows as the number of samples in \code{measurementsTrain} and 2 to 3 columns if it is a
#' \code{\link{matrix}} or a \code{\link{DataFrame}}, or a character vector of length 2 to 3 containing the
#' column names in \code{measurementsTrain} if it is a \code{\link{DataFrame}} or the
#' column name in \code{colData(measurementsTrain)} if \code{measurementsTrain} is a
#' \code{\link{MultiAssayExperiment}}. If a vector of column names, those columns will be
#' removed before training.
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
#' @importFrom survival coxph
#' @rdname coxphRanking
#' @export
setGeneric("coxphRanking", function(measurementsTrain, ...) standardGeneric("coxphRanking"))

#' @rdname coxphRanking
#' @export
setMethod("coxphRanking", "matrix", function(measurementsTrain, survivalTrain, ...) # Matrix of numeric measurements.
{
  coxphRanking(S4Vectors::DataFrame(measurementsTrain, check.names = FALSE), survivalTrain, ...)
})

#' @rdname coxphRanking
#' @export
setMethod("coxphRanking", "DataFrame", function(measurementsTrain, survivalTrain, verbose = 3) # Clinical data or one of the other inputs, transformed.
{
  splitDataset <- .splitDataAndOutcomes(measurementsTrain, survivalTrain)
  measurementsTrain <- splitDataset[["measurements"]]
  survivalTrain <- splitDataset[["outcomes"]]

  pValues <- apply(measurementsTrain, 2, function(featureColumn){
    fit <- survival::coxph(survivalTrain ~ featureColumn)
    s <- summary(fit)
    s$waldtest["pvalue"]
  })
  
  if(!is.null(S4Vectors::mcols(measurementsTrain)))
    S4Vectors::mcols(measurementsTrain)[order(pValues), ]
  else
    colnames(measurementsTrain)[order(pValues)]
})

# One or more omics data sets, possibly with clinical data.
#' @rdname coxphRanking
#' @export
setMethod("coxphRanking", "MultiAssayExperiment", function(measurementsTrain, targets = names(measurementsTrain), survivalTrain, ...)
{
  tablesAndSurvival <- .MAEtoWideTable(measurementsTrain, targets, survivalTrain)
  measurementsTrain <- tablesAndSurvival[["dataTable"]]
  survivalTrain <- tablesAndSurvival[["outcomes"]]
  
  if(ncol(measurementsTrain) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    coxphRanking(measurementsTrain, survivalTrain, ...)
})