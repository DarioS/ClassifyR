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
  \S4method{runTest}{ExpressionSet}(expression, doFirst = c("transform", "selection"),
                     training, testing, transformParams, selectionParams, trainParams,
                     predictParams, verbose = 1)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    the training data. For a matrix, the rows are features, and the columns
                    are samples.}
  \item{classes}{A vector of class labels.}
  \item{doFirst}{Character. Either "selection" or "transform". Specifies whether the first
                 step of the process should be selection or transformation.}
  \item{training}{A vector which specifies the training samples.}
  \item{testing}{A vector which specifies the test samples.}    
  \item{transformParams}{A container of class \code{\link{TransformParams}} which contains
                         the function to do transformation and parameters it will use.}
  \item{selectionParams}{A container of class \code{\link{SelectionParams}} which contains
                         the function to do feature selection and parameters it will use.}       
  \item{trainParams}{A container of class \code{\link{TrainParams}} describing the
                     classifier to use for training.}
  \item{predictParams}{A container of class \code{\link{PredictParams}} describing how
                       prediction is to be done.}
  \item{...}{Unused variables from the \code{\link{matrix}} method passed to the
             \code{\link{ExpressionSet}} method.}
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 A higher number will produce more messages.}
}
\details{
  This function only performs one classification and prediction. See \code{\link{runTests}}
  for a driver function that does cross validation and uses this function.
}
\value{
  A list with three elements. The first element contains the indices of genes that were selected
  by the feature selection step. The second element contains the indices of the samples that
  were in the test set. The third element contains a vector of the classes predicted by the
  classifer.
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
    runTest(TCGA_eset, doFirst = "selection", training = (1:ncol(TCGA_eset)) \%\% 2 == 0,
            testing = (1:ncol(TCGA_eset)) \%\% 2 != 0)
  }
}
