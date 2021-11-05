\name{TrainParams}
\docType{class}
\alias{TrainParams}
\alias{TrainParams-class}
\alias{TrainParams,missing-method}
\alias{TrainParams,function-method}
\alias{show,TrainParams-method}

\title{Parameters for Classifier Training}
\description{
  Collects and checks necessary parameters required for classifier training. The 
  empty constructor is provided for convenience.
}
\section{Constructor}{
  \describe{
    \item{}{
      \code{TrainParams()}
            Creates a default TrainParams object. The classifier function is \code{dlda} for
            Diagonal LDA. Users should create an appropriate \code{TrainParams} object for the
            characteristics of their data, once they are familiar with this software.
  }
    \item{}{
      \preformatted{TrainParams(classifier, characteristics = DataFrame(), intermediate = character(0),
                                getFeatures = NULL, ...)}
            Creates a TrainParams object which stores the function which will do the
            classifier building and parameters that the function will use.
  \describe{
            \item{\code{classifier}}{A function which will construct a classifier, and also
            possibly make the predictions. The first argument must be a \code{\link{DataFrame}}
            object. The second argument must be a vector of classes. If the function also makes predictions
            and the value of the \code{predictor} setting of \code{PredictParams} is therefore \code{NULL},
            the third argument must be a \code{DataFrame} of test data. The function must also accept a
            parameter named \code{verbose}. The function's return value can be either a trained classifier
            if the function only does training or a vector or data frame of class predictions if
            it also does prediction with the test set samples.}
            \item{\code{characteristics}}{A \code{\link{DataFrame}} describing the characteristics
            of the classifier used. First column must be named \code{"charateristic"} and second column must
            be named \code{"value"}. If using wrapper functions for classifiers
            in this package, a classifier name will automatically be generated and therefore it is
            not necessary to specify it.}
            \item{\code{intermediate}}{Character vector. Names of any variables created in prior stages by
            \code{\link{runTest}} that need to be passed to \code{classifier}.}
            \item{\code{getFeatures}}{A function may be specified that extracts the selected features from
            the trained model. This is relevant if using a classifier that does feature selection within
            training (e.g. random forest). The function must return a list of two vectors. The first vector
            contains the ranked features (or empty if the training algorithm doesn't produce rankings) and
            the second vector contains the selected features.}
            \item{\code{...}}{Other named parameters which will be used by the classifier.}
	         }
    }
  }
}

\section{Summary}{
  A method which summarises the parameter settings.
  \code{trainParams} is a \code{TrainParams} object.

  \describe{
    \item{}{
      \code{show(trainParams)}: Prints a short summary of what \code{trainParams} contains.
    }}
}

\author{Dario Strbenac}
\examples{
#if(require(sparsediscrim))
  trainParams <- TrainParams(DLDAtrainInterface)
}