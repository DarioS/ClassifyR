\name{previousTrained}
\alias{previousTrained}
\alias{previousTrained,ClassifyResult-method}
\title{Automated Usage of Previously Created Classifiers}
\description{Uses the trained classifier of the same cross-validation iteration of a previous classification for the current classification task.
}
\usage{
  \S4method{previousTrained}{ClassifyResult}(classifyResult, .iteration, verbose = 3)
}
\arguments{
  \item{classifyResult}{A \code{\link{ClassifyResult}} object which stores the models fitted previously.}
  \item{.iteration}{Do not specify this variable. It is set by \code{\link{runTests}} if this function
                    is being repeatedly called by \code{runTests}.}
  \item{verbose}{Default: 3. A number between 0 and 3 for the amount of progress messages to give.
    This function only prints progress messages if the value is 3.}                    
}
\value{
  A trained classifier from a previously completed classification task.
}
\author{Dario Strbenac}

\examples{
  #if(require(sparsediscrim))
  #{
    # Genes 76 to 100 have differential expression.
    genesMatrix <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
    genesMatrix <- cbind(genesMatrix, sapply(1:25, function(sample)
                                 c(rnorm(75, 9, 2), rnorm(25, 14, 2))))
    colnames(genesMatrix) <- paste("Sample", 1:50)
    rownames(genesMatrix) <- paste("Gene", 1:100)                                 
    classes <- factor(rep(c("Poor", "Good"), each = 25))
    resubstitute <- ResubstituteParams(nFeatures = seq(10, 100, 10),
                                       performanceType = "error")
    result <- runTests(genesMatrix, classes,
                       permutations = 2, fold = 2,
                       params = list(SelectParams(), TrainParams(), PredictParams()))
    models(result)
                       
    # Genes 50 to 74 have differential expression in new data set.
    newDataset <- sapply(1:25, function(sample) c(rnorm(100, 9, 2)))
    newDataset <- cbind(newDataset, rbind(sapply(1:25, function(sample) rnorm(49, 9, 2)),
                                          sapply(1:25, function(sample) rnorm(25, 14, 2)),
                                          sapply(1:25, function(sample) rnorm(26, 9, 2))))
                                          
    rownames(newDataset) <- rownames(genesMatrix)
    colnames(newDataset) <- colnames(genesMatrix)
                                           
    newerResult <- runTests(newDataset, classes,
                            permutations = 2, fold = 2,
                            params = list(SelectParams(previousSelection,
                                                       intermediate = ".iteration",
                                                       classifyResult = result),
                                          TrainParams(previousTrained,
                                                      intermediate = ".iteration",
                                                      classifyResult = result),
                                          PredictParams()))
                                     
    models(newerResult)
  #}  
}
