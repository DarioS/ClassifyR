\name{SelectParams}
\docType{class}
\alias{SelectParams}
\alias{SelectParams-class}
\alias{SelectParams,missing-method}
\alias{SelectParams,functionOrList-method}

\title{Parameters for Feature Selection}
\description{
  Collects and checks necessary parameters required for feature selection. Either one
  function is specified or a list of functions to perform ensemble feature selection. The 
  empty constructor is provided for convenience.
}
\section{Constructor}{
  \describe{
    \item{}{
      \code{SelectParams()}
      Creates a default SelectParams object. This uses either an ordinary t-test or ANOVA
      (depending on the number of classes) and tries the top 10 to top 100 features in
      increments of 10, and picks the number of features with the best resubstitution
      balanced error rate. Users should create an appropriate \code{SelectParams}
      object for the characteristics of their data, once they are familiar with this software.
  }
    \item{}{
      \preformatted{SelectParams(featureSelection, characteristics = DataFrame(),
      minPresence = 1, intermediate = character(0), subsetToSelections = TRUE, ...)}
      Creates a \code{SelectParams} object which stores the function which will do the
      selection and parameters that the function will use.
	\describe{
            \item{\code{featureSelection}}{Either a function which will do the selection
            or a list of such functions. For a particular function, the first argument must
            be an \code{\link{DataFrame}} object. The function's return value must be
            a \code{\link{SelectResult}} object.}
            \item{\code{characteristics}}{A \code{\link{DataFrame}} describing the characteristics
            of feature selection to be done. First column must be named \code{"charateristic"} and
            second column must be named \code{"value"}. If using wrapper functions for feature selection
            in this package, the feature selection name will automatically be generated and therefore it is
            not necessary to specify it.}
            \item{\code{minPresence}}{If a list of functions was provided, how many of those
            must a feature have been selected by to be used in classification. 1 is equivalent
            to a set union and a number the same length as \code{featureSelection}
            is equivalent to set intersection.}
            \item{\code{intermediate}}{Character vector. Names of any variables created in prior stages by
            \code{\link{runTest}} that need to be passed to a feature selection function.}
            \item{\code{subsetToSelections}}{Whether to subset the data table(s), after
                                             feature selection has been done.}
	          \item{\code{...}}{Other named parameters which will be used by the selection function.
            If \code{featureSelection} was a list of functions, this must be a list of lists, as
            long as \code{featureSelection}.}
	         }
    }
  }
}

\author{Dario Strbenac}
\examples{
  #if(require(sparsediscrim))
  #{
    SelectParams(differentMeansSelection,
                 trainParams = TrainParams(), predictParams = PredictParams(),
                 resubstituteParams = ResubstituteParams())
    
    # Ensemble feature selection.
    SelectParams(list(differentMeansSelection, pairsDifferencesSelection),
                 trainParams = TrainParams(), predictParams = PredictParams(),
                 resubstituteParams = ResubstituteParams())
  #}
}