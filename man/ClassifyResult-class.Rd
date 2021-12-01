\name{ClassifyResult}
\docType{class}
\alias{ClassifyResult}
\alias{ClassifyResult-class}
\alias{ClassifyResult,DataFrame,character,characterOrDataFrame-method}
\alias{show,ClassifyResult-method}
\alias{sampleNames}
\alias{sampleNames,ClassifyResult-method}
\alias{featureNames}
\alias{featureNames,ClassifyResult-method}
\alias{predictions}
\alias{predictions,ClassifyResult-method}
\alias{actualClasses}
\alias{actualClasses,ClassifyResult-method}
\alias{features}
\alias{features,ClassifyResult-method}
\alias{models}
\alias{models,ClassifyResult-method}
\alias{performance}
\alias{performance,ClassifyResult-method}
\alias{tunedParameters}
\alias{tunedParameters,ClassifyResult-method}
\alias{totalPredictions}
\alias{totalPredictions,ClassifyResult-method}

\title{Container for Storing Classification Results}
\description{
  Contains a list of models, table of actual sample classes and predicted classes,
  the identifiers of features selected for each fold of each permutation or
  each hold-out classification, and performance metrics such as error rates.
  This class is not intended to be created by the user. It is created by \code{\link{runTest}} or \code{\link{runTests}}.
}
\section{Constructor}{
  \describe{
    \item{}{
      \preformatted{ClassifyResult(characteristics, originalNames, originalFeatures,
           rankedFeatures, chosenFeatures, predictions, actualClasses, models,
           validation, tune = NULL)}
            }
    }
  \describe{
            \item{\code{characteristics}}{A \code{\link{DataFrame}} describing the characteristics
            of classification done. First column must be named \code{"charateristic"} and second column must
            be named \code{"value"}. If using wrapper functions for feature selection and classifiers
            in this package, the function names will automatically be
            generated and therefore it is not necessary to specify them.}
            \item{\code{originalNames}}{All sample names.}
            \item{\code{originalFeatures}}{All feature names. Character vector or \code{\link{DataFrame}}
                                           with one row for each feature if the data set is a
                                          \code{\link{MultiAssayExperiment}}.}
            \item{\code{rankedFeatures}}{All features, from most to least important. Character vector
                                         or \code{\link{DataFrame}} if data set is a
                                         \code{MultiAssayExperiment}.}
            \item{\code{chosenFeatures}}{Features selected at each fold. Character vector or DataFrame if
                                         data set is a \code{MultiAssayExperiment}.}
            \item{\code{predictions}}{A \code{\link{list}} of \code{\link{data.frame}}
                                      containing information about samples, their actual class and
                                      predicted class and/or class score and information about the
                                      cross-validation fold in which the prediction was made.}            
            \item{\code{actualClasses}}{Factor of class of each sample.}                                      
            \item{\code{models}}{All of the models fitted to the training data.}                                    
            \item{\code{validation}}{List with first element being the name of the validation scheme,
            and other elements providing details about the scheme.}
            \item{\code{tune}}{A description of the tuning parameters, and the value chosen of
                               each parameter.}
           }
}

\section{Summary}{
  A method which summarises the results is available.
  \code{result} is a \code{ClassifyResult} object.

  \describe{
    \item{}{
      \code{show(result)}: Prints a short summary of what \code{result} contains.
    }}
  \describe{
    \item{}{
      \code{totalPredictions(ClassifyResult)}: Calculates the sum of the number of predictions.
    }}  
}

\section{Accessors}{
  \code{result} is a \code{ClassifyResult} object.

  \describe{
  \item{\code{sampleNames(result)}}{Returns a vector of sample names present in the data set.}}
  \describe{
  \item{\code{featureNames(result)}}{Returns a vector of features present in the data set.}}
  \describe{
  \item{\code{predictions(result)}}{Returns a \code{list} of \code{data.frame}.
     Each data.frame contains columns \code{sample}, \code{predicted}, and \code{actual}. For
     hold-out validation, only one data.frame is returned of all of the concatenated
     predictions.}}
  \describe{
  \item{\code{actualClasses(result)}}{Returns a \code{factor} class labels, one for
    each sample.}}
  \describe{
  \item{\code{features(result)}}{A \code{list} of the features selected for each training.}}
  \describe{
  \item{\code{models(result)}}{A \code{list} of the models fitted for each training.}}
  \describe{
  \item{\code{performance(result)}}{Returns a \code{list} of performance measures. This is
    empty until \code{calcCVperformance} has been used.}}
  \describe{
  \item{\code{tunedParameters(result)}}{Returns a \code{list} of tuned parameter values.
        If cross-validation is used, this list will be large, as it stores chosen values
        for every iteration.}}    
  \describe{
  \item{\code{sampleNames(result)}}{Returns a \code{\link{character}} vector of sample names.}}
}

\author{Dario Strbenac}
\examples{
  #if(require(sparsediscrim))
  #{
    data(asthma)
    
    LOOCVparams <- CrossValParams("Leave-k-Out", leave = 1)
    modellingParams <- ModellingParams()
    classified <-
    runTests(measurements, classes, LOOCVparams, modellingParams,
             DataFrame(characteristic = c("dataset", "classification"),
                      value = c("Asthma", "Different Means"))
             )
    class(classified)
  #}
}
