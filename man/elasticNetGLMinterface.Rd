\name{elasticNetGLMinterface}
\alias{elasticNetGLMinterface}
\alias{elasticNetGLMtrainInterface}
\alias{elasticNetGLMpredictInterface}
\alias{elasticNetGLMtrainInterface,matrix-method}
\alias{elasticNetGLMtrainInterface,DataFrame-method}
\alias{elasticNetGLMtrainInterface,MultiAssayExperiment-method}
\alias{elasticNetGLMpredictInterface,multnet,matrix-method}
\alias{elasticNetGLMpredictInterface,multnet,DataFrame-method}
\alias{elasticNetGLMpredictInterface,multnet,MultiAssayExperiment-method}
\title{
  An Interface for glmnet Package's glmnet Function
}
\description{
  An elastic net GLM classifier uses a penalty which is a combination of a lasso penalty and
  a ridge penalty, scaled by a lambda value, to fit a sparse linear model to the data.
}
\usage{
  \S4method{elasticNetGLMtrainInterface}{matrix}(measurements, classes, ...)
  \S4method{elasticNetGLMtrainInterface}{DataFrame}(measurements, classes, lambda = NULL,
                             ..., verbose = 3)
  \S4method{elasticNetGLMtrainInterface}{MultiAssayExperiment}(measurements, targets = names(measurements), ...)
  \S4method{elasticNetGLMpredictInterface}{multnet,matrix}(model, test, ...)
  \S4method{elasticNetGLMpredictInterface}{multnet,DataFrame}(model, test, classes = NULL, lambda, ..., returnType = c("class", "score", "both"), verbose = 3)
  \S4method{elasticNetGLMpredictInterface}{multnet,MultiAssayExperiment}(model, test, targets = names(test), ...)                                   
}
\arguments{
  \item{measurements}{Either a \code{\link{matrix}}, \code{\link{DataFrame}} or
                     \code{\link{MultiAssayExperiment}} containing the training data.
                     For a \code{matrix}, the rows are features, and the columns are samples.
                     If of type \code{\link{DataFrame}}, the data set is subset to only those
                     features of type \code{integer}.}
  \item{classes}{Either a vector of class labels of class \code{\link{factor}} of the same length
                 as the number of samples in \code{measurements} or if the measurements are
                 of class \code{DataFrame} a character vector of length 1 containing the
                 column name in \code{measurement} is also permitted. Not used if \code{measurements}
                 is a \code{MultiAssayExperiment} object.}
  \item{lambda}{The lambda value passed directly to \code{\link[glmnet]{glmnet}} if the training
  function is used or passed as \code{s} to \code{\link[glmnet]{predict.glmnet}} if the prediction
  function is used.}
  \item{test}{An object of the same class as \code{measurements} with no samples in common with
              \code{measurements} and the same number of features as it.}                
  \item{targets}{If \code{measurements} is a \code{MultiAssayExperiment}, the names of the
                 data tables to be used. \code{"clinical"} is also a valid value and specifies that
                 integer variables from the clinical data table will be used.}
  \item{...}{Variables not used by the \code{matrix} nor the \code{MultiAssayExperiment} method which
             are passed into and used by the \code{DataFrame} method (e.g. \code{verbose}) or, for the
             training function, options that are used by the \code{glmnet} function. For the testing 
             function, this variable simply contains any parameters passed from the classification
             framework to it which aren't used by glmnet's \code{predict} fuction.}
  \item{model}{A trained elastic net GLM, as created by the \code{glmnet} function.}
  \item{returnType}{Default: \code{"class"}. Either \code{"class"}, \code{"score"} or \code{"both"}.
                    Sets the return value from the prediction to either a vector of class labels,
                    matrix of scores for each class, or both labels and scores in a \code{data.frame}.}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
                 This function only prints progress messages if the value is 3.}             
}
\details{
  If \code{measurements} is an object of class \code{MultiAssayExperiment}, the factor of sample
  classes must be stored in the DataFrame accessible by the \code{colData} function with
  column name \code{"class"}.
  
  The value of the \code{"family"} parameter is fixed to \code{"multinomial"} so that
  classification with more than 2 classes is possible. During classifier training, if more than
  one lambda value is considered by specifying a vector of them as input, then the chosen value
  is determined based on classifier resubstitution performance. Leaving the default value of NULL,
  however, causes the \code{\link[glmnet]{glmnet}} function to determine a set of values to evaluate
  and then a single lambda value is chosen based on the fraction of (null) deviance explained by the
  fitted model.
}
\value{
  For \code{elasticNetGLMtrainInterface}, an object of type \code{glmnet}. For
  \code{elasticNetGLMpredictInterface}, either a factor vector of predicted classes, a
  matrix of scores for each class, or a table of both the class labels and class scores,
  depending on the setting of \code{returnType}.
}
\author{Dario Strbenac}
\examples{
  if(require(glmnet))
  {
    # Genes 76 to 100 have differential expression.
    genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
                                      c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    colnames(genesMatrix) <- paste("Sample", 1:ncol(genesMatrix))
    rownames(genesMatrix) <- paste("Gene", 1:nrow(genesMatrix))
    
    resubstituteParams <- ResubstituteParams(nFeatures = seq(10, 100, 10),
                                             performanceType = "balanced error",
                                             better = "lower")
                                             
    # alpha is a user-specified tuning parameter.
    # lambda is automatically tuned, based on glmnet defaults, if not user-specified.                 
    trainParams <- TrainParams(elasticNetGLMtrainInterface,
                               tuneParams = list(alpha = c(0, 0.5, 1)),
                               tuneOptimise = c("balanced error", "lower"))
    predictParams <- PredictParams(elasticNetGLMpredictInterface)                           
    classified <- runTests(genesMatrix, classes, datasetName = "Example",
                           classificationName = "Differential Expression",
                           validation = "leaveOut", leave = 1,
                           params = list(trainParams, predictParams))
                           
    classified <- calcCVperformance(classified, "balanced error")
    head(tunedParameters(classified))
    performance(classified)
  }
}
