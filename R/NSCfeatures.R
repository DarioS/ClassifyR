#' Interface for \code{pamr.listgenes} Function from \code{pamr} CRAN Package
#' 
#' Extracts the threshold for the minimum training error and then extracts the
#' corresponding gene IDs of the genes that were not eliminated by the
#' thresold.
#' 
#' When used within ClassifyR cross-validation, the trained model, measurements
#' and classes will automatically be passed to this function in each iteration.
#' 
#' @aliases NSCfeatures NSCfeatures,pamrtrained-method
#' @param trained The output of \code{\link{NSCtrainInterface}}, which is
#' identical to the output of \code{\link[pamr]{pamr.listgenes}}.
#' @param measurements A \code{\link{DataFrame}} containing the training data.
#' @param classes A vector of class labels of class \code{\link{factor}} of the
#' same length as the number of samples in \code{measurements}.
#' @return A list with the first element being empty (no feature ranking is
#' provided) and second element being the selected features.
#' @author Dario Strbenac
#' @seealso \code{\link[pamr]{pamr.listgenes}} for the function that is
#' interfaced to.
#' @examples
#' 
#'   if(require(pamr))
#'   {
#'     # Genes 76 to 100 have differential expression.
#'     genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(100, 9, 1)))
#'     genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn)
#'                                  c(rnorm(75, 9, 1), rnorm(25, 14, 1))))
#'     rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))                                 
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     
#'     trained <- NSCtrainInterface(genesMatrix, classes)
#'     # ClassifyR framework internally uses DataFrames for measurements storage.
#'     selected <- NSCfeatures(trained, DataFrame(t(genesMatrix), check.names = FALSE), classes)
#'     selected[[2]]                                                       
#'   }
#' 
setGeneric("NSCfeatures", function(trained, measurements, classes)
standardGeneric("NSCfeatures"))

setMethod("NSCfeatures", "pamrtrained",
          function(trained, measurements, classes)
{
  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")
  
  minError <- min(trained[["errors"]])
  threshold <- trained[["threshold"]][max(which(trained[["errors"]] == minError))]
  params <- c(list(trained), list(list(x = t(as.matrix(measurements)), y = classes, geneid = 1:ncol(measurements))), threshold)
  chosen <- as.numeric(do.call(pamr::pamr.listgenes, params)[, 1])
  
  if(is.null(S4Vectors::mcols(measurements)))
    chosen <- colnames(measurements)[chosen]
  else
    chosen <- S4Vectors::mcols(measurements)[chosen, ]
  
  list(NULL, chosen)
})