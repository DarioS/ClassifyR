\name{elasticNetFeatures}
\alias{elasticNetFeatures}
\alias{elasticNetFeatures,multnet-method}
\title{Extract Vectors of Ranked and Selected Features From an Elastic Net GLM Object}
\description{Provides a ranking of features based on the magnitude of fitted GLM coefficients. Also provides the selected features which are those with a non-zero coefficient.
}
\usage{
  \S4method{elasticNetFeatures}{multnet}(model)
}
\arguments{
  \item{model}{A fitted multinomial GLM which was created by \code{\link[glmnet]{glmnet}}.}
}
\value{
  An \code{list} object. The first element is a vector or data frame of ranked features, the second is a vector or data frame of selected features.
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
                                               performanceType = "balanced error")
                                             
      # alpha is a user-specified tuning parameter.
      # lambda is automatically tuned, based on glmnet defaults, if not user-specified.                 
      trainParams <- TrainParams(elasticNetGLMtrainInterface, nlambda = 500)
      predictParams <- PredictParams(elasticNetGLMpredictInterface)                           
      classified <- runTests(genesMatrix, classes, validation = "fold",
                             params = list(trainParams, predictParams))
                                        
      elasticNetFeatures(models(classified)[[1]])
    }
}
