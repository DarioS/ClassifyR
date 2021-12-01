\name{StageParams}
\docType{class}
\alias{StageParams}
\alias{StageParams-class}

\title{StageParams Virtual Class}
\description{
  A class for any one of \code{\link{TransformParams}}, \code{\link{SelectParams}}, \code{\link{TrainParams}} or \code{\link{PredictParams}}. Allows a method to dispatch
  on any of the paramter objects specifying any stage of cross-validation.
}

\author{Dario Strbenac}