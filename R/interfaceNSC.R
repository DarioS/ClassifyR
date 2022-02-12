################################################################################
#
# Train interface
#
################################################################################

#' Interface for \code{pamr.train} Function from \code{pamr} CRAN Package
#' 
#' Restructures variables from ClassifyR framework to be compatible with
#' \code{\link[pamr]{pamr.train}} definition.
#' 
#' This function is an interface between the ClassifyR framework and
#' \code{\link[pamr]{pamr.train}}.
#' 
#' @aliases NSCtrainInterface NSCtrainInterface,matrix-method
#' NSCtrainInterface,DataFrame-method
#' NSCtrainInterface,MultiAssayExperiment-method
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
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
#' @param measurementsTest An object of the same class as \code{measurementsTrain} with no
#' samples in common with \code{measurementsTrain} and the same number of features
#' as it.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method or extra arguments passed to \code{\link[pamr]{pamr.train}}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return A list with elements as described in \code{\link[pamr]{pamr.train}}.
#' @author Dario Strbenac
#' @seealso \code{\link[pamr]{pamr.train}} for the function that was interfaced
#' to.
#' @examples
#' 
#'   if(require(pamr))
#'   {
#'     # Samples in one class with differential expression to other class.
#'     genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(100, 9, 1)))
#'     genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn)
#'                                  c(rnorm(75, 9, 1), rnorm(25, 14, 1))))
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     
#'     NSCtrainInterface(genesMatrix, classes)
#'   }
#' 
#' @export
setGeneric("NSCtrainInterface", function(measurementsTrain, ...)
standardGeneric("NSCtrainInterface"))

setMethod("NSCtrainInterface", "matrix", function(measurementsTrain, classes, ...)
{
  NSCtrainInterface(DataFrame(measurementsTrain, check.names = FALSE), classes, ...)
})

setMethod("NSCtrainInterface", "DataFrame", # Sample information data or one of the other inputs, transformed.
          function(measurementsTrain, classesTrain, ..., verbose = 3)
{
  splitDataset <- .splitDataAndOutcomes(measurementsTrain, classesTrain)
  measurementsTrain <- splitDataset[["measurements"]]
  classesTrain <- splitDataset[["outcomes"]]

  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")

  # The features are expected to be in the rows and the samples are in the columns.
  trainedModel <- pamr::pamr.train(list(x = t(as.matrix(measurementsTrain)), y = classesTrain), ...)
  
  if(verbose == 3)
    message("Nearest shrunken centroid training completed.")
  
  trainedModel  
})

setMethod("NSCtrainInterface", "MultiAssayExperiment",
          function(measurementsTrain, targets = names(measurementsTrain), classesTrain, ...)
{ 
  tablesAndClasses <- .MAEtoWideTable(measurementsTrain, targets, classesTrain)
  measurementsTrain <- tablesAndClasses[["dataTable"]]
  classesTrain <- tablesAndClasses[["outcomes"]]
  
  if(ncol(measurementsTrain) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    NSCtrainInterface(measurementsTrain, classesTrain, ...)
})

################################################################################
#
# Predict interface
#
################################################################################


#' Interface for \code{pamr.predict} Function from \code{pamr} CRAN Package
#' 
#' Restructures variables from ClassifyR framework to be compatible with
#' \code{\link[pamr]{pamr.predict}} definition.
#' 
#' This function is an interface between the ClassifyR framework and
#' \code{\link[pamr]{pamr.predict}}.  It selects the highest threshold that
#' gives the minimum error rate in the training data.
#' 
#' @aliases NSCpredictInterface NSCpredictInterface,pamrtrained,matrix-method
#' NSCpredictInterface,pamrtrained,DataFrame-method
#' NSCpredictInterface,pamrtrained,MultiAssayExperiment-method
#' @param model An object of class \code{pamrtrained}.
#' @param measurementsTest An object of the same class as \code{measurementsTrain} with no
#' samples in common with \code{measurementsTrain} used in the training stage and
#' the same number of features as it.  Also, if a \code{DataFrame}, the
  #' \code{classesTrain} column must be absent.
#' @param classesColumnTest Either NULL or a character vector of length 1, specifying the
#' column name to remove from the test set.
#' @param targets If \code{measurementsTest} is a \code{MultiAssayExperiment}, the names of
#' the data tables to be used. \code{"sampleInfo"} is also a valid value and
#' specifies that numeric variables from the sample information table will be used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method or optional settings that are passed to
#' \code{\link[pamr]{pamr.predict}}.
#' @param returnType Default: \code{"both"}. Either \code{"class"},
#' \code{"score"} or \code{"both"}.  Sets the return value from the prediction
#' to either a vector of class labels, score for a sample belonging to the
#' second class, as determined by the factor levels, or both labels and scores
#' in a \code{data.frame}.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return Either a factor vector of predicted classes, a matrix of scores for
#' each class, or a table of both the class labels and class scores, depending
#' on the setting of \code{returnType}.
#' @author Dario Strbenac
#' @seealso \code{\link[pamr]{pamr.predict}} for the function that was
#' interfaced to.
#' @examples
#' 
#'   if(require(pamr))
#'   {
#'     # Samples in one class with differential expression to other class.
#'     genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(100, 9, 1)))
#'     genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn)
#'                                  c(rnorm(75, 9, 1), rnorm(25, 14, 1))))
#'     classes <- factor(rep(c("Poor", "Good"), each = 25))
#'     
#'     fit <- NSCtrainInterface(genesMatrix[, c(1:20, 26:45)], classes[c(1:20, 26:45)])
#'     NSCpredictInterface(fit, genesMatrix[, c(21:25, 46:50)])
#'   }
#' 
#' @export
setGeneric("NSCpredictInterface", function(model, measurementsTest, ...)
  standardGeneric("NSCpredictInterface"))

setMethod("NSCpredictInterface", c("pamrtrained", "matrix"), function(model, measurementsTest, ...)
{
  NSCpredictInterface(model, DataFrame(test, check.names = FALSE), ...)
})

setMethod("NSCpredictInterface", c("pamrtrained", "DataFrame"), function(model, measurementsTest, classesColumnTest = NULL, ..., returnType = c("both", "class", "score"), verbose = 3)
{
  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")
  returnType <- match.arg(returnType)
  
  if(!is.null(classesColumnTest)) # Remove the column, since pamr uses positional matching of features.
  {
    splitDataset <- .splitDataAndOutcomes(measurementsTest, classesColumnTest) 
    measurementsTest <- splitDataset[["measurements"]] # Without classes column.
  }
  
  minError <- min(model[["errors"]])
  threshold <- model[["threshold"]][max(which(model[["errors"]] == minError))]
  
  measurementsTest <- t(as.matrix(measurementsTest))   
  classPredictions <- pamr::pamr.predict(model, measurementsTest, threshold, ...)
  classScores <- pamr::pamr.predict(model, measurementsTest, threshold, type = "posterior", ...)[, levels(model[["y"]])]
  if(!is.matrix(classScores)) # Only one sample was predicted and pamr isn't consistent with return types.
    classScores <- t(classScores)
  
  if(verbose == 3)
    message("Nearest shrunken centroid predictions made.")
  
  switch(returnType, class = classPredictions, # Factor vector.
         score = classScores, # Numeric matrix.
         both = data.frame(class = classPredictions, classScores, check.names = FALSE))
})

setMethod("NSCpredictInterface", c("pamrtrained", "MultiAssayExperiment"), function(model, measurementsTest, targets = names(measurementsTest), ...)
{
  measurementsTest <- .MAEtoWideTable(measurementsTest, targets)[["dataTable"]] # Remove any classes, if present.
  
  if(ncol(measurementsTest) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    NSCpredictInterface(model, measurementsTest, ...)
})


################################################################################
#
# Get selected features
#
################################################################################

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
#' @param model The output of \code{\link{NSCtrainInterface}}, which is
#' identical to the output of \code{\link[pamr]{pamr.listgenes}}.
#' @param measurementsTrain A \code{\link{DataFrame}} containing the training data.
#' @param classesTrain A vector of class labels of class \code{\link{factor}} of the
#' same length as the number of samples in \code{measurementsTrain}.
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
#'     model <- NSCtrainInterface(genesMatrix, classes)
#'     # ClassifyR framework internally uses DataFrames for measurements storage.
#'     selected <- NSCfeatures(model, DataFrame(t(genesMatrix), check.names = FALSE), classes)
#'     selected[[2]]                                                       
#'   }
#' 
setGeneric("NSCfeatures", function(model, measurementsTrain, classesTrain)
  standardGeneric("NSCfeatures"))

setMethod("NSCfeatures", "pamrtrained",
          function(model, measurementsTrain, classesTrain)
          {
            if(!requireNamespace("pamr", quietly = TRUE))
              stop("The package 'pamr' could not be found. Please install it.")
            
            minError <- min(model[["errors"]])
            threshold <- model[["threshold"]][max(which(model[["errors"]] == minError))]
            params <- c(list(model), list(list(x = t(as.matrix(measurementsTrain)), y = measurementsTrain, geneid = 1:ncol(measurementsTrain))), threshold)
            chosen <- as.numeric(do.call(pamr::pamr.listgenes, params)[, 1])
            
            if(is.null(S4Vectors::mcols(measurementsTrain)))
              chosen <- colnames(measurementsTrain)[chosen]
            else
              chosen <- S4Vectors::mcols(measurementsTrain)[chosen, ]
            
            list(NULL, chosen)
          })