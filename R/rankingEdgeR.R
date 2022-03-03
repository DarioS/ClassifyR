#' Feature Ranking Based on Differential Expression for Count Data
#' 
#' Performs a differential expression analysis between classes and ranks the
#' features based on test statistics. The data may have overdispersion and this
#' is modelled.
#' 
#' The differential expression analysis follows the standard
#' \code{\link[edgeR]{edgeR}} steps of estimating library size normalisation
#' factors, calculating dispersion, in this case robustly, and then fitting a
#' generalised linear model followed by a likelihood ratio test.
#' 
#' Data tables which consist entirely of non-numeric data cannot be analysed.
#' 
#' @aliases edgeRranking edgeRranking,matrix-method
#' edgeRranking,DataFrame-method edgeRranking,MultiAssayExperiment-method
#' @param countsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
#' \code{\link{MultiAssayExperiment}} containing the unnormalised counts.
#' For a \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns
#' are features, unlike the convention used in edgeR.
#' @param classesTrain A vector of class labels of class \code{\link{factor}} of the
#' same length as the number of samples in \code{countsTrain} if it is a
#' \code{\link{matrix}} or a \code{\link{DataFrame}} or a character vector of
#' length 1 containing the column name in \code{countsTrain} if it is a
#' \code{\link{DataFrame}} or the column name in \code{colData(counts)} if
#' \code{countsTrain} is a \code{\link{MultiAssayExperiment}}. If a column name,
#' that column will be removed before training.
#' @param targets If \code{countsTrain} is a \code{MultiAssayExperiment}, the
#' names of the data tables of counts to be used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method.
#' @param normFactorsOptions A named \code{list} of any options to be passed to
#' \code{\link[edgeR]{calcNormFactors}}.
#' @param dispOptions A named \code{list} of any options to be passed to
#' \code{\link[edgeR]{estimateDisp}}.
#' @param fitOptions A named \code{list} of any options to be passed to
#' \code{\link[edgeR]{glmFit}}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return A vector or data frame (if \code{MultiAssayExperiment} input) of
#' features, from the most promising features in the first position to the
#' least promising feature in the last position.
#' @author Dario Strbenac
#' @references edgeR: a Bioconductor package for differential expression
#' analysis of digital gene expression data, Mark D. Robinson, Davis McCarthy,
#' and Gordon Smyth, 2010, \emph{Bioinformatics}, Volume 26 Issue 1,
#' \url{https://academic.oup.com/bioinformatics/article/26/1/139/182458}.
#' @examples
#' 
#'   if(require(parathyroidSE) && require(PoiClaClu))
#'   {
#'     data(parathyroidGenesSE)
#'     expression <- assays(parathyroidGenesSE)[[1]] # Genes in rows, samples in columns.
#'     sampleNames <- paste("Sample", 1:ncol(parathyroidGenesSE))
#'     colnames(expression) <- sampleNames
#'     DPN <- which(colData(parathyroidGenesSE)[, "treatment"] == "DPN")
#'     control <- which(colData(parathyroidGenesSE)[, "treatment"] == "Control")
#'     expression <- expression[, c(control, DPN)]
#'     classes <- factor(rep(c("Contol", "DPN"), c(length(control), length(DPN))))
#'     expression <- expression[rowSums(expression > 1000) > 8, ] # Make small data set.
#'     
#'     # ClassifyR is using the convention of samples in rows and features in columns.
#'     ranked <- edgeRranking(t(expression), classes)
#'                                         
#'     head(ranked)
#'     plotFeatureClasses(t(expression), classes, "ENSG00000044574",
#'                        dotBinWidth = 500, xAxisLabel = "Unnormalised Counts")
#'   }
#' 
#' @export
setGeneric("edgeRranking", function(countsTrain, ...) standardGeneric("edgeRranking"))

setMethod("edgeRranking", "matrix", function(countsTrain, classesTrain, ...) # Matrix of integer counts. 
{
  edgeRranking(DataFrame(countsTrain, check.names = FALSE), classesTrain, ...)
})

# DataFrame of counts, likely created by runTests or runTest.
setMethod("edgeRranking", "DataFrame", function(countsTrain, classesTrain, normFactorsOptions = NULL, dispOptions = NULL, fitOptions = NULL, verbose = 3)
{
  if(verbose == 3)
    message("Doing edgeR LRT feature ranking")
  
  # DGEList stores features as rows and samples as columns.          
  countsList <- edgeR::DGEList(t(as.matrix(countsTrain)), group = classesTrain)
  paramList <- list(countsList)
  if(!is.null(normFactorsOptions))
    paramList <- append(paramList, normFactorsOptions)
  if(verbose == 3)
    message("Calculating scaling factors.")
  countsList <- do.call(edgeR::calcNormFactors, paramList)
  paramList <- list(countsList, model.matrix(~ classesTrain))
  if(!is.null(dispOptions))
    paramList <- append(paramList, dispOptions)
  if(verbose == 3)
    message("Estimating dispersion.")
  countsList <- do.call(edgeR::estimateDisp, paramList)
  paramList <- list(countsList, model.matrix(~ classesTrain))
  if(!is.null(fitOptions))
    paramList <- append(paramList, fitOptions)
  if(verbose == 3)
    message("Fitting linear model.")
  fit <- do.call(edgeR::glmFit, paramList)
  test <- edgeR::glmLRT(fit, coef = 2:length(levels(classesTrain)))[["table"]]
  
  if(!is.null(S4Vectors::mcols(countsTrain)))
    S4Vectors::mcols(countsTrain)[order(test[, "PValue"]), ]
  else
    colnames(countsTrain)[order(test[, "PValue"])]
})

# One or more omics data sets, possibly with sample information data.
setMethod("edgeRranking", "MultiAssayExperiment", function(countsTrain, targets = NULL, ...)
{
  if(!requireNamespace("edgeR", quietly = TRUE))
    stop("The package 'edgeR' could not be found. Please install it.")
  if(is.null(targets))
    stop("'targets' must be specified but was not.")
  if(length(setdiff(targets, names(countsTrain))))
    stop("Some values of 'targets' are not names of 'counts' but all must be.")            

  tablesAndClasses <- .MAEtoWideTable(countsTrain, targets, "integer")
  countsTable <- tablesAndClasses[["dataTable"]]
  classesTrain <- tablesAndClasses[["outcomes"]]
  edgeRranking(countsTable, classesTrain, ...)
})