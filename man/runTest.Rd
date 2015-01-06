\name{runTest}
\alias{runTest}
\alias{runTest,matrix-method}
\alias{runTest,ExpressionSet-method}
\title{Perform a Single Classification}
\description{For a dataset of features and samples, the classification process is run.
It consists of data transformation, feature selection, training and testing.
       }
\usage{
  \S4method{runTest}{matrix}(expression, classes, ...)
  \S4method{runTest}{ExpressionSet}(expression, datasetName, classificationName,
                     training, testing, params = list(SelectionParams(), TrainParams(), PredictParams()),
                     verbose = 1)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    the training data. For a matrix, the rows are features, and the columns
                    are samples.}
  \item{classes}{A vector of class labels.}
  \item{datasetName}{A name associated with the dataset used.}
  \item{classificationName}{A name associated with the classification.}  
  \item{training}{A vector which specifies the training samples.}
  \item{testing}{A vector which specifies the test samples.}    
  \item{params}{A \code{\link{list}} of objects of class of \code{\link{TransformParams}},
               \code{\link{SelectionParams}}, \code{\link{TrainParams}}, or \code{\link{PredictParams}}.
               The order they are in the list determines the order in which the stages
               of classification are done in.}
  \item{...}{Unused variables from the \code{\link{matrix}} method passed to the
             \code{\link{ExpressionSet}} method.}
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 A higher number will produce more messages.}
}
\details{
  This function only performs one classification and prediction. See \code{\link{runTests}}
  for a driver function that does cross validation and uses this function. \code{datasetName}
  and \code{classificationName} need to be provided.
}
\value{
  A list with five elements. The first element contained all of the features, ranked from most
  important to least important.The second element contains the indices of genes that were selected
  by the feature selection step. The third element contains the indices of the samples that
  were in the test set. The fourth element contains a vector of the classes predicted by the
  classifer. The fifth element contains the value of any tuning parameters tried and chosen.
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
    runTest(TCGA_eset, "TCGA Cancer", "Differential Expression",
            training = (1:ncol(TCGA_eset)) \%\% 2 == 0,
            testing = (1:ncol(TCGA_eset)) \%\% 2 != 0)
  }
}
