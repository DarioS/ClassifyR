\name{ClassifyResult}
\docType{class}
\alias{ClassifyResult}
\alias{ClassifyResult-class}
\alias{ClassifyResult,character,character,character,character-method}
\alias{show,ClassifyResult-method}
\alias{predictions}
\alias{predictions,ClassifyResult-method}
\alias{actualClasses}
\alias{actualClasses,ClassifyResult-method}
\alias{features}
\alias{features,ClassifyResult-method}
\alias{performance}
\alias{performance,ClassifyResult-method}
\alias{tunedParameters}
\alias{tunedParameters,ClassifyResult-method}
\alias{totalPredictions}
\alias{totalPredictions,ClassifyResult-method}

\title{Container for Storing Classification Results}
\description{
  Contains a table of actual sample classes and predicted classes, the indices of
  features selected for each fold of each bootstrap resampling or each hold-out
  classification, and error rates. This class is not intended to be created by
  the user, but could be used in another package. It is created by \code{\link{runTests}}.
}
\section{Constructor}{
  \describe{
    \item{}{
      \code{ClassifyResult(datasetName, classificationName, originalNames, originalFeatures, rankedFeatures,
      chosenFeatures, predictions, actualClasses, validation, tune = list(NULL))}}
    }
  \describe{
            \item{\code{datasetName}}{A name associated with the dataset used.}
            \item{\code{classificationName}}{A name associated with the classification.}
            \item{\code{originalNames}}{Sample names.}
            \item{\code{originalFeatures}}{Feature names.}
            \item{\code{rankedFeatures}}{Indices or names of all features, from most to least
                                         important.}
            \item{\code{chosenFeatures}}{Indices or names of features selected at each fold.}
            \item{\code{predictions}}{A \code{\link{list}} of \code{\link{data.frame}}
                                      containing information about samples, their actual class and
                                      predicted class.}            
            \item{\code{actualClasses}}{Factor of class of each sample.}
            \item{\code{validation}}{List with first elment being name of the validation scheme,
            and other elements providing details about scehme.}
            \item{\code{tune}}{A description of the tuning parameters, and the value chosen of
                               each parameter.}
           }
}

\section{Summary}{
  A method which summarises the results is available.
  \code{result} is a \code{ClassifyResult} object.

  \describe{
    \item{}{
      \code{show(result)}{Prints a short summary of what \code{result} contains.}
    }}
  \describe{
    \item{}{
      \code{totalPredictions(ClassifyResult)}{Calculates the sum of the number of predictions.}
    }}  
}

\section{Accessors}{
  \code{result} is a \code{ClassifyResult} object.

  \describe{
  \item{\code{predictions(result)}}{Returns a \code{\link{list}} of \code{\link{data.frame}}.
     Each data.frame contains columns \code{sample}, \code{predicted}, and \code{actual}. For
     hold-out validation, only one data.frame is returned of all of the concatenated
     predictions.}}
  \describe{
  \item{\code{actualClasses(result)}}{Returns a \code{\link{factor}} class labels, one for
    each sample.}}
  \describe{
  \item{\code{features(result)}}{Returns a \code{\link{list}} of \code{\link{data.frame}}.
    Each data.frame contains columns \code{sample}, \code{predicted}, and \code{actual}.}}
  \describe{
  \item{\code{performance(result)}}{Returns a \code{\link{list}} of performance measures. This is
    empty until \code{\link{calcPerformance}} has been used.}}
  \describe{
  \item{\code{tunedParameters(result)}}{Returns a \code{\link{list}} of tuned parameter values. If cross-validation is used, this list will be large, as it stores chosen values for every validation.}}    
  \describe{
  \item{\code{names(result)}}{Returns a \code{\link{character}} vector of sample names.}}
}

\author{Dario Strbenac}
\examples{
  if(require(curatedOvarianData) && require(sparsediscrim))
  {
    data(TCGA_eset)
    badOutcome <- which(pData(TCGA_eset)[, "vital_status"] == "deceased" & pData(TCGA_eset)[, "days_to_death"] <= 365)
    goodOutcome <- which(pData(TCGA_eset)[, "vital_status"] == "living" & pData(TCGA_eset)[, "days_to_death"] >= 365 * 5)
    TCGA_eset <- TCGA_eset[, c(badOutcome, goodOutcome)]
    classes <- factor(rep(c("Poor", "Good"), c(length(badOutcome), length(goodOutcome))))
    pData(TCGA_eset)[, "class"] <- classes
    results <- runTests(TCGA_eset, "Ovarian Cancer", "Differential Expression", resamples = 2, folds = 2)
    show(results)
    predictions(results)
    actualClasses(results)
  }
}
