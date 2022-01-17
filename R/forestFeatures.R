#' Extract Vectors of Ranked and Selected Features From a Random Forest Object
#'
#' Provides a ranking of features based on the total decrease in node
#' impurities from splitting on the variable, averaged over all trees. Also
#' provides the selected features which are those that were used in at least
#' one tree of the forest.
#'
#'
#' @aliases forestFeatures forestFeatures,randomForest-method
#' @param forest A trained random forest which was created by
#' \code{\link{randomForest}}.
#' @return An \code{list} object. The first element is a vector or data frame
#' of features, ranked from best to worst using the Gini index. The second
#' element is a vector or data frame of features used in at least one tree.
#' @author Dario Strbenac
#' @examples
#'
#'       if(require(randomForest))
#'       {
#'         genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
#'         genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
#'                                           c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
#'         classes <- factor(rep(c("Poor", "Good"), each = 25))
#'         colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix))
#'         rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
#'         trainingSamples <- c(1:20, 26:45)
#'         testingSamples <- c(21:25, 46:50)
#'
#'         trained <- randomForestTrainInterface(genesMatrix[, trainingSamples],
#'                                               classes[trainingSamples], ntree = 10)
#'
#'         forestFeatures(trained)
#'       }
#'
#' @importFrom randomForest importance varUsed
#' @export
setGeneric("forestFeatures", function(forest, ...)
           standardGeneric("forestFeatures"))

setMethod("forestFeatures", "randomForest",
          function(forest)
{
  inputFeatures <- rownames(randomForest::importance(forest))
  rankedFeatures <- inputFeatures[order(randomForest::importance(forest), decreasing = TRUE)]
  selectedFeatures <- inputFeatures[randomForest::varUsed(forest) > 0]
  selectedFeatures <- selectedFeatures[na.omit(match(rankedFeatures, selectedFeatures))]

  # Colon is a reserved symbol for separating data name and feature name, which is necessary for identifiability of MultiAssayExperiment features. It is not permitted in feature names.
  if(grepl(':', selectedFeatures[1]) == TRUE) # Convert to data.frame.
  {
    selectedFeatures <- do.call(rbind, strsplit(selectedFeatures, ':'))
    rankedFeatures <- do.call(rbind, strsplit(rankedFeatures, ':'))
    colnames(selectedFeatures) <- colnames(rankedFeatures) <- c("dataset", "feature")
  }
  list(rankedFeatures, selectedFeatures)
})