\name{easyHardFeatures}
\alias{easyHardFeatures}
\alias{easyHardFeatures,EasyHardClassifier-method}
\title{
  Extract Chosen Features from an EasyHardClassifier Object
}
\description{
  The features are described by a data frame. One column is named "dataset" and the other is named "feature". This provides identifiability in case when multiple types of data have features with the same name.
}
\usage{
  \S4method{easyHardFeatures}{EasyHardClassifier}(easyHardClassifier)
}
\arguments{
  \item{easyHardClassifier}{An \code{\link{EasyHardClassifier}} object}.
}
\details{
  If any of the features are from the colun data of the input \code{\link{MultiAssayExperiment}}, the dataset value will be "clinical".
}
\value{
  To be consistent with other functions for extracting features from a trained model, a list of length two. The first element is for feature rankings, which is not meaningful for an easy-hard classifier, so it is \code{NULL}. The second element is the selected features.
}
\author{Dario Strbenac}
\examples{
  genesMatrix <- matrix(c(rnorm(90, 9, 1),
                          9.5, 9.4, 5.2, 5.3, 5.4, 9.4, 9.6, 9.9, 9.1, 9.8),
		        ncol = 10, byrow = TRUE)
  colnames(genesMatrix) <- paste("Sample", 1:10)
  rownames(genesMatrix) <- paste("Gene", 1:10)
  genders <- factor(c("Male", "Male", "Female", "Female", "Female",
                      "Female", "Female", "Female", "Female", "Female"))

  # Scenario: Male gender can predict the hard-to-classify Sample 1.
  clinical <- DataFrame(age = c(31, 34, 32, 39, 33, 38, 34, 37, 35, 36),
                        gender = genders,
                        class = factor(rep(c("Poor", "Good"), each = 5)),
		        row.names = colnames(genesMatrix))
  dataset <- MultiAssayExperiment(ExperimentList(RNA = genesMatrix), clinical)
  trained <- easyHardClassifierTrain(dataset, easyClassifierParams = list(minCardinality = 2, minPurity = 0.9),
                                     hardClassifierParams = list(SelectParams(featureSelection = differentMeansSelection,
                                                                 selectionName = "Difference in Means",
                                                                 resubstituteParams = ResubstituteParams(1:10, "balanced error", "lower")), TrainParams(), PredictParams()))

  easyHardFeatures(trained)
}
