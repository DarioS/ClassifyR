\name{DMDselection}
\alias{DMDselection}
\alias{DMDselection,matrix-method}
\alias{DMDselection,ExpressionSet-method}
\title{Selection of Differential Distributions with Kullback Leibler Distance}
\description{Ranks features by largest Differences in Means/Medians and Deviations and chooses
the features which have best resubstitution performance.
       }
\usage{
  \S4method{DMDselection}{matrix}(expression, classes, ...)
  \S4method{DMDselection}{ExpressionSet}(expression, trainParams,
                                         predictParams, resubstituteParams, ..., verbose = 3)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    the training data. For a matrix, the rows are features, and the columns
                    are samples.}
  \item{classes}{A vector of class labels.}
  \item{trainParams}{A container of class \code{\link{TrainParams}} describing the
                     classifier to use for training.}
  \item{predictParams}{A container of class \code{\link{PredictParams}} describing how
                       prediction is to be done.}
  \item{resubstituteParams}{An object of class \code{\link{ResubstituteParams}}
                            describing the performance measure to consider and the numbers of
                            top features to try for resubstitution classification.}                        
  \item{...}{Variables passed to \code{\link{getLocationsAndScales}}.}                       
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}
}
\details{
  DMD is defined as \eqn{|location_1 - location_2| + |scale_1 -  scale_2|}{|location1 - location2| + |scale1 -  scale2|}.
  
  The subscripts denote the group which the parameter is calculated for.
}
\value{
  A list of length 2. The first element has the features ranked from most important to least important.
  The second element has the features that were selected to be used for classification.
}
\author{Dario Strbenac}

\examples{
  if(require(sparsediscrim))
  {
    # First 25 samples are mixtures of two normals. Last 25 samples are one normal.
    genesMatrix <- sapply(1:25, function(geneColumn) c(rnorm(50, 5, 1), rnorm(50, 15, 1)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(geneColumn) rnorm(100, 9, 3)))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    DMDselection(genesMatrix, classes,
                 trainParams = TrainParams(), predictParams = PredictParams(),
                 resubstituteParams = ResubstituteParams(nFeatures = seq(10, 100, 10), performanceType = "balanced", better = "lower"))
  }
}