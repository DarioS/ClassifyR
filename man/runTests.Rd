\name{runTests}
\alias{runTests}
\alias{runTests,matrix-method}
\alias{runTests,ExpressionSet-method}
\title{Reproducibly Do Resampling or Leave Out and Cross Validation}
\description{Enables doing classification schemes such as 100 resamples 5-fold cross validation
             or leave one out cross validaion. Processing in parallel is possible by leveraging
             the package \code{\link{BiocParallel}}.
}
\usage{
  \S4method{runTests}{matrix}(expression, classes, ...)
  \S4method{runTests}{ExpressionSet}(expression, datasetName, classificationName,
                     validation = c("bootstrap", "leaveOut"), bootMode = c("fold", "split"),
                     resamples = 100, percent = 25, folds = 5, leave = 2, seed, parallelParams = bpparam(),
                     params = list(SelectionParams(), TrainParams(), PredictParams()),
                     verbose = 1)
}
\arguments{
  \item{expression}{Either a \code{\link{matrix}} or \code{\link{ExpressionSet}} containing
                    the training data. For a matrix, the rows are features, and the columns
                    are samples.}
  \item{classes}{A vector the same length as the number of columns of expression data specifying
                 the class that the samples belong to.}            
  \item{datasetName}{A name associated with the dataset used.}
  \item{classificationName}{A name associated with the classification.}
  \item{validation}{"bootstrap" for repeated resampling or "leaveOut" for leaving all combinations
                    of k samples as test samples.}
  \item{bootMode}{Character. Either "fold" or "split". If "fold", then the samples are split into
                 folds and in each iteration one is used as the test set. If "split", the samples
                 are split into two groups. One is used as the training set, the other is the test set.}                    
  \item{resamples}{Relevant when repeated resampling is used. The number of times to do
                   sampling with replacement.}
  \item{percent}{Used when bootstrap resampling with split method is chosen. The percentage of
                 samples to be in the test set.}                   
  \item{folds}{Relevant when repeated resampling is used with fold mode. The number of folds
               to break each resampling into. Each fold is used once as the test set.}
  \item{leave}{Relevant when leave k out validation is used. The number of samples to leave
               for testing.}
  \item{seed}{The random number generator used for repeated resampling will use this seed,
              if it is provided. Allows reproducibility of repeated usage on the same input data.}
  \item{parallelParams}{An object of class \code{\link{MulticoreParam}} or \code{\link{SnowParam}}.}               
  \item{params}{A \code{\link{list}} of objects of class of \code{\link{TransformParams}},
               \code{\link{SelectionParams}}, \code{\link{TrainParams}}, or \code{\link{PredictParams}}.
               The order they are in the list determines the order in which the stages
               of classification are done in.}
  \item{...}{Unused variables from the \code{\link{matrix}} method passed to the
             \code{\link{ExpressionSet}} method.}
  \item{verbose}{A number between 0 and 3 for the amount of progress messages to give.
                 A higher number will produce more messages.}
}
\value{
  If the predictor function made a single prediction, then an object of class
  \code{\link{ClassifyResult}}. If the predictor function made a set of predictions, then
  a list of such objects.
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
    runTests(TCGA_eset, "Ovarian Cancer", "Differential Expression", resamples = 2, fold = 2)
  }
}
