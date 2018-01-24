\name{logisticRegressionInterface}
\alias{logisticRegressionTrainInterface}
\alias{logisticRegressionTrainInterface,matrix-method}
\alias{logisticRegressionTrainInterface,DataFrame-method}
\alias{logisticRegressionTrainInterface,MultiAssayExperiment-method}
\alias{logisticRegressionPredictInterface}
\alias{logisticRegressionPredictInterface,mlogit,matrix-method}
\alias{logisticRegressionPredictInterface,mlogit,DataFrame-method}
\alias{logisticRegressionPredictInterface,mlogit,MultiAssayExperiment-method}
\title{
  An Interface for mlogit Package's mlogit Function
}
\description{
  \code{logisticRegressionTrainInterface} generates a multinomial logistic regression model trained on some
  training data and \code{logisticRegressionPredictInterface} makes class predictions for samples in 
  the test dataset.
}
\usage{
  \S4method{logisticRegressionTrainInterface}{matrix}(measurements, classes, ...)
  \S4method{logisticRegressionTrainInterface}{DataFrame}(measurements, classes, ...)
  \S4method{logisticRegressionTrainInterface}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
  \S4method{logisticRegressionPredictInterface}{matrix}(model, test, ...)
  \S4method{logisticRegressionPredictInterface}{DataFrame}(model, test, ...)
  \S4method{logisticRegressionPredictInterface}{MultiAssayExperiment}(model, test, targets = names(measurements), ...)  
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
                     \code{\link{MultiAssayExperiment}} containing the training data.
                     For a \code{matrix}, the rows are features, and the columns are samples.}
  \item{classes}{Either a vector of class labels of class \code{\link{factor}} of the same length
                 as the number of samples in \code{measurements} or if the measurements are
                 of class \code{DataFrame} a character vector of length 1 containing the
                 column name in \code{measurement} is also permitted. Not used if \code{measurements}
                 is a \code{MultiAssayExperiment} object.}
  \item{test}{An object of the same class as \code{measurements} with no samples in common with
              \code{measurements} and the same number of features as it.}                 
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 integer variables from the clinical data table will be used.}
  \item{...}{Unused variables by the three top-level methods passed to the internal method
             which does the training or prediction (i.e. \code{verbose}).}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}  
}
\details{
  If \code{measurements} is an object of class \code{MultiAssayExperiment}, the factor of
  sample classes must be stored in the DataFrame accessible by the \code{colData} function
  with column name \code{"class"}.
}
\value{
  For \code{logisticRegressionTrainInterface}, a fitted multinomial logistic regression model.
  For \code{logisticRegressionPredictInterface}, a result of type \code{list} with elements \code{class},
  \code{scores} and \code{posterior}, as created by \code{sparsediscrim}'s \code{predict} method.
}
\author{Dario Strbenac}
\examples{
  # 3 cancer types, melanoma, breast cancer, prostate cancer.
  # 30 samples, 10 samples have each cancer type.
  
  cancerTypes <- c("Melanoma", "Breast Cancer", "Prostate Cancer")
  clinical <- DataFrame(Age = rpois(30, c(rep(40,5), rep(50, 5), rep(35, 20))),
                        Gender = factor(c(rep(c("Male", "Female"), c(5, 5)),
                                          rep(c("Male", "Female"), c(1, 9)),
                                          rep("Male", 10))),
                        Cancer = factor(rep(cancerTypes, each = 10)),
                        row.names = paste("Patient", 1:30))
                        
  geneExpression <- matrix(c(rnorm(30, 15, 1), # Housekeeping gene.
                             c(rnorm(10, 12, 1), rnorm(20, 7, 1)), # High in melanoma.
                             c(rnorm(10, 8, 1), rnorm(10, 9, 1), rnorm(10, 7, 2)), # High in breast cancer.
                             c(rnorm(20, 10, 1), rnorm(10, 6, 2)) # Low in prostate cancer.
                             ), ncol = 30, byrow = TRUE,
                           dimnames = list(row = c("GAPDH", "BRAF", "FST", "GSTP1"),
                                           column = paste("Patient", 1:30))
                          )
  datasetContainer <- MultiAssayExperiment(list(RNA = geneExpression), clinical)
}
