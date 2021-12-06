\name{ModellingParams}
\docType{class}
\alias{ModellingParams}
\alias{ModellingParams-class}

\title{Parameters for Data Modelling Specification}
\description{
  Collects and checks necessary parameters required for data modelling. Apart from data transfomation that needs to be
  done within cross-validation (e.g. \code{\link{subtractFromLocation}}), feature selection,
  model training and prediction, this container also stores a setting for class imbalance rebalancing.

}
\section{Constructor}{
  \describe{
    \item{}{
      \code{ModellingParams()}
      Creates a default \code{ModellingParams} object suitable for classification of homoscedastic numeric data.
      This uses either an ordinary t-test or ANOVA (depending on the number of classes) and tries the top 10 to
      top 100 features in increments of 10 during feature selection and diagonal LDA as the classifier. Also,
      the default class rebalancing scheme is downsampling to the smallest class. Users should create an
      appropriate \code{ModellingParams} object for the characteristics of their data.
  }
    \item{}{
      \preformatted{ModellingParams(balancing = c("downsample", "upsample", "none"), transformParams = NULL,
                   selectParams = SelectParams(), trainParams = TrainParams(), predictParams = PredictParams())}
      Creates a \code{ModellingParams} object based on user-specified settings with default values of
      different classification stages geared at homoscedastic numeric data.
	\describe{
            \item{\code{balancing}}{Default: "downsample". A character value specifying what kind of class balancing to do, if any.}
            \item{\code{transformParams}}{Parameters used for feature transformation specified by a \code{\link{TransformParams}} instance.
            Optional, can be \code{NULL}}.
            \item{\code{selectParams}}{Parameters used during feature selection specified by a \code{\link{SelectParams}} instance.
            By default, parameters for selection based on differences in means of numeric data. Optional, can be \code{NULL}.}
            \item{\code{trainParams}}{Parameters for model training specified by a \code{\link{TrainParams}} instance.
            By default, uses diagonal LDA.}
            \item{\code{predictParams}}{Parameters for model training specified by a \code{\link{PredictParams}} instance.
            By default, uses diagonal LDA.}
	         }
    }
  }
}

\author{Dario Strbenac}
\examples{
  #if(require(sparsediscrim))
  #{
     ModellingParams() # Default is differences in means selection and DLDA.
     ModellingParams(selectParams = NULL, # No feature selection before training.
                     trainParams = TrainParams(randomForestTrainInterface),
                     predictParams = PredictParams(randomForestPredictInterface))
  #}
}