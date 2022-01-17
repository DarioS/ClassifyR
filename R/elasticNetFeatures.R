#' Extract Vectors of Ranked and Selected Features From an Elastic Net GLM
#' Object
#' 
#' Provides a ranking of features based on the magnitude of fitted GLM
#' coefficients. Also provides the selected features which are those with a
#' non-zero coefficient.
#' 
#' 
#' @aliases elasticNetFeatures elasticNetFeatures,multnet-method
#' @param model A fitted multinomial GLM which was created by
#' \code{\link[glmnet]{glmnet}}.
#' @return An \code{list} object. The first element is a vector or data frame
#' of ranked features, the second is a vector or data frame of selected
#' features.
#' @author Dario Strbenac
#' @examples
#' 
#'     if(require(glmnet))
#'     {
#'       # Genes 76 to 100 have differential expression.
#'       genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
#'       genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
#'                                         c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
#'       classes <- factor(rep(c("Poor", "Good"), each = 25))
#'       colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix))
#'       rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
#'                                              
#'       # alpha is a user-specified tuning parameter.
#'       # lambda is automatically tuned, based on glmnet defaults, if not user-specified.                 
#'       CVparams <- CrossValParams("k-Fold")
#'       
#'       trainParams <- TrainParams(elasticNetGLMtrainInterface, nlambda = 500)
#'       predictParams <- PredictParams(elasticNetGLMpredictInterface)
#'       modParams <- ModellingParams(selectParams = NULL, trainParams = trainParams,
#'                                    predictParams = predictParams)
#'                                    
#'       classified <- runTests(genesMatrix, classes, CVparams, modParams)
#'                                         
#'       elasticNetFeatures(models(classified)[[1]])
#'     }
#' 
#' @export
setGeneric("elasticNetFeatures", function(model, ...)
           standardGeneric("elasticNetFeatures"))

setMethod("elasticNetFeatures", "multnet",
          function(model)
{
  inputFeatures <- rownames(model[["beta"]][[1]])            
  # Floating point numbers test for equality.
  whichCoefficientColumn <- which(abs(model[["lambda"]] - attr(model, "tune")[["lambda"]]) < 0.00001)[1]
  coefficientsUsed <- sapply(model[["beta"]], function(classCoefficients) classCoefficients[, whichCoefficientColumn])
  featureScores <- rowSums(abs(coefficientsUsed))
  rankedFeatures <- inputFeatures[order(featureScores, decreasing = TRUE)]
  selectedFeatures <- inputFeatures[featureScores != 0]
  
  # Colon is a reserved symbol for separating data name and feature name, which is necessary for identifiability of MultiAssayExperiment features. It is not permitted in feature names.
  if(grepl(':', selectedFeatures[1]) == TRUE) # Convert to data.frame.
  {
    selectedFeatures <- do.call(rbind, strsplit(selectedFeatures, ':'))
    rankedFeatures <- do.call(rbind, strsplit(rankedFeatures, ':'))
    colnames(selectedFeatures) <- colnames(rankedFeatures) <- c("dataset", "feature")
  }
  list(unique(rankedFeatures), selectedFeatures)
})