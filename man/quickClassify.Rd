\name{quickClassify}
\alias{quickClassify}
\alias{quickClassify,matrix-method}
\alias{quickClassify,DataFrame-method}
\alias{quickClassify,MultiAssayExperiment-method}
\title{Quick and Easy Cross-validated Classification}
\description{This convenience function allows specification of cross-validation and classification
by a couple of pre-defined keywords and all parameter setting is automatically taken care of. This avoids the need to create any S4 parameter objects.
}
\usage{
  \S4method{quickClassify}{matrix}(measurements, classes, ...)
  \S4method{quickClassify}{DataFrame}(measurements, classes,
crossVal = c("100 permutations, 5 folds", "100 permutations, 25\% test", "5-fold", "leave-1-out", "leave-2-out"),
                   seed = 12345, cores = 8,
classifier = c("DLDA", "FishersLDA", "NSC", "naiveBayes", "elasticNet", "SVM", "randomForest", "kTSP"),
                  ...)
  \S4method{quickClassify}{MultiAssayExperiment}(measurements, targets, ...)
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
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 numeric variables from the clinical data table will be used.}
  \item{crossVal}{Default: "100 permutations, 5 folds". A short name describing the kind of cross-validation to do. Self-explanatory.}
  \item{seed}{Default: 12345. Random number seed to use. Sample splitting and classification on parallel threads will be reproducible.}
  \item{cores}{Default: 8. Number of CPU cores to use. Each cross-validation iteration is done on one core.}
  \item{classifier}{Default: "DLDA". A short name describing the kind of classification to do. The options represent:
  \itemize{\item \code{"DLDA"}: Feature selection based on differences in means ranking and diagonal LDA classifier.
  \item \code{"FishersLDA"}: Transformation by subtraction of values from the mean value of all samples in the training set followed by selection by Levene's test statistic ranking based on differential variability and classification by Fisher's LDA.
  \item \code{"NSC"}: Nearest shrunken centroid classifier. No explicit feature selection before classifier.
  \item \code{"naiveBayes"}: Feature selection based on differences in median and Qn followed by naive Bayes classifier using weighted voting based on the difference between class densities. See \code{\link{DMDranking}} for more details. Good for differential dsitribution classification.
    \item \code{"elasticNet"}: Elastic Net GLM classifier with \code{"multinomial"} as the family.
    \item \code{"SVM"}: Polynomial kernel support vector machine. Does parameter tuning of poly nomial degree between 2 and 8 and the cost penality between 0.00001 to 100000 in power increments of 1.
    \item \code{"randomForest"}: Random forest with the number of variables tried at each split eing 40\% of the total number of variables in the data set.
    \item \code{"kTSP"}: Feature selection based on pairs differences ranking and classification based on weighted voting of the selected feature pairs.
  }
  }
  \item{...}{Either variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method or variables not used by \code{quickClassify}
             but used by runTests (i.e. \code{characteristics} and \code{verbose}.}
}
\details{
  The relevant cross-validation and modelling parameters are created and \code{\link{runTests} is executed}. Parallel processing is reproducible because the seed is used for splitting samples into training and test sets, but also for the parallel loop which runs each iteration of cross-validation.
  
  For shortcuts which do feature selection, the top 10, 20, ..., 100 ranked features are evaluated for resubstitution error rate unless the data set has less than 100 features, in which case, all possible top sets of \emph{p} features are evaluated (i.e. 1, 2, 3, ..., p - 1, p) before one of them is chosen.
}
\value{
  An object of class \code{\link{ClassifyResult}}.
}
\author{Dario Strbenac}

\examples{
  data(asthma) # Restrict to 20 samples for speed.
  measurements <- measurements[, 1:20]
  classes <- classes[1:20]
  result <- quickClassify(measurements, classes, crossVal = "leave-1-out",
                          classifier = "randomForest", seed = 99)
}
