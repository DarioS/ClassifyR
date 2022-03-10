#' Automated Usage of Previously Created Classifiers
#' 
#' Uses the trained classifier of the same cross-validation iteration of a
#' previous classification for the current classification task.
#' 
#' 
#' @aliases previousTrained previousTrained,ClassifyResult-method
#' @param classifyResult A \code{\link{ClassifyResult}} object which stores the
#' models fitted previously.
#' @param .iteration Do not specify this variable. It is set by
#' \code{\link{runTests}} if this function is being repeatedly called by
#' \code{runTests}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return A trained classifier from a previously completed classification
#' task.
#' @author Dario Strbenac
#' @examples
#' 
#'   #if(require(sparsediscrim))
#'   #{
#'     # Genes 76 to 100 have differential expression.
#'     genesMatrix <- sapply(1:100, function(sample) rnorm(25, 9, 0.3))
#'     genesMatrix <- rbind(genesMatrix, t(sapply(1:25, function(sample)
#'                                       c(rnorm(75, 9, 0.3), rnorm(25, 14, 0.3)))))
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     rownames(genesMatrix) <- paste("Sample", 1:50)
#'     colnames(genesMatrix) <- paste("Gene", 1:100)                                                           
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#' 
#'     CVparams <- CrossValParams(permutations = 2, folds = 2)
#'     result <- runTests(genesMatrix, classes, CVparams, ModellingParams())
#'     models(result)
#'                        
#'     # Genes 50 to 74 have differential expression in new data set.
#'     newDataset <- sapply(1:100, function(sample) c(rnorm(25, 9, 0.3)))
#'     newDataset <- rbind(newDataset, cbind(sapply(1:49, function(sample) rnorm(25, 9, 0.3)),
#'                                           sapply(1:25, function(sample) rnorm(25, 14, 0.3)),
#'                                           sapply(1:26, function(sample) rnorm(25, 9, 0.3))))
#'                                           
#'     rownames(newDataset) <- rownames(genesMatrix)
#'     colnames(newDataset) <- colnames(genesMatrix)
#' 
#'     selPars <- SelectParams(previousSelection, intermediate = ".iteration", classifyResult = result)
#'     trPars <- TrainParams(previousTrained, intermediate = ".iteration", classifyResult = result)
#'     previousParams <- ModellingParams(selectParams = selPars, trainParams = trPars)                           
#'     newerResult <- runTests(newDataset, classes, CVparams, previousParams)
#'     models(newerResult)
#'   #}  
#'
#' @usage NULL
#' @export
setGeneric("previousTrained", function(classifyResult, ...)
standardGeneric("previousTrained"))

#' @export
setMethod("previousTrained", "ClassifyResult", 
          function(classifyResult, .iteration, verbose = 3)
{
  if(verbose == 3)
    message("Using existing classifier for classification.")
  
  models(classifyResult)[[.iteration]]
})