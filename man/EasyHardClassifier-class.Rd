\name{EasyHardClassifier}
\docType{class}
\alias{EasyHardClassifier}
\alias{EasyHardClassifier-class}
\alias{EasyHardClassifier,listOrNULL,listOrCharacterOrNULL,character-method}
\alias{show,EasyHardClassifier-method}

\title{Container for a Pair of Trained Classifiers}
\description{
  Stores two classifiers and the ID of the data set each one of them was trained on. Not intended for end-user.
}
\section{Constructor}{
  \describe{
    \item{}{
      \preformatted{EasyHardClassifier(easyClassifier, hardClassifier, datasetIDs)}
            Creates a EasyHardClassifier object which stores the details of the two underlying classification models.
  \describe{
            \item{\code{easyClassifier}}{A classifier trained on the easy-to-collect data set by \link{easyHardClassifierTrain}. Represented as a list of rules, each of which is also a list.}
            \item{\code{hardClassifier}}{A list of two with names "selected" and "model". The "selected" element should contain
            the chosen features and the "model" element is for the trained model on the hard-to-collect data set or simply a character vector of length 1 containing a class name, if the prediction of the samples left over from the easy classifier all or all except one belonged to a particular class.}
            \item{\code{datasetIDs}}{A vector of length 2 with names "easy" and "hard" containing the data set IDs from the original \code{\link{MultiAssayExperiment}} input data object for the easy data set and the hard data set.}
	         }
    }
  }
}

\section{Summary}{
  A method which shows the easy and hard classifiers.
  \code{easyHard} is a \code{EasyHardClassifier} object.

  \describe{
    \item{}{
      \code{show(EasyHardClassifier)}: Prints a short summary of what \code{easyHard} contains.
    }}
}

\author{Dario Strbenac}
\examples{
  predictiveRules <- list(list(feature = "age", relation = "<", value = 18, predict = "Risk"))
  hardClassifier <- DLDAtrainInterface(matrix(rnorm(400), ncol = 20),
                        classes = factor(sample(c("Safe", "Risk"), 20, replace = TRUE)))
  EasyHardClassifier(predictiveRules, list(selected = LETTERS[1:5], model = hardClassifier),
                     setNames(c("clinical", "RNA-seq"), c("easy", "hard")))
}