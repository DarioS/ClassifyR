\name{PredictParamsOrNULL}
\docType{class}
\alias{PredictParamsOrNULL}
\alias{PredictParamsOrNULL-class}

\title{Union of A PredictParams Object and NULL}
\description{
  Allows a slot to be either a PredictParams class object or empty. No constructor. In other words,
  the training function specified also makes predictions with the test set.
}

\author{Dario Strbenac}

\examples{
  ModellingParams(trainParams = TrainParams(kNNinterface, k = 5), predictParams = NULL)
}