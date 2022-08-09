#' Ranking of Differential Variability with coxph Statistic
#' 
#' Ranks all features from largest coxph statistic to smallest.
#' 
#' The calculation of the test statistic is performed by the
#' coxph function from the survival package.
#' 
#' Data tables which consist entirely of non-numeric data cannot be ranked.
#' 
#' @aliases coxphRanking coxphRanking,matrix-method
#' coxphRanking,DataFrame-method coxphRanking,MultiAssayExperiment-method
#' @param measurementsTrain Either a \code{\link{matrix}}, \code{\link{DataFrame}}
#' or \code{\link{MultiAssayExperiment}} containing the training data. For a
#' \code{matrix} or \code{\link{DataFrame}}, the rows are samples, and the columns are features.
#' @param survivalTrain A tabular data type of survival information of the
#' same number of rows as the number of samples in \code{measurementsTrain} and 2 to 3 columns if it is a
#' \code{\link{matrix}} or a \code{\link{DataFrame}}, or a character vector of length 2 to 3 containing the
#' column names in \code{measurementsTrain} if it is a \code{\link{DataFrame}} or the
#' column name in \code{colData(measurementsTrain)} if \code{measurementsTrain} is a
#' \code{\link{MultiAssayExperiment}}. If a vector of column names, those columns will be
#' removed before training.
#' @param targets If \code{measurements} is a \code{MultiAssayExperiment}, the
#' names of the data tables to be used. \code{"clinical"} is also a valid value
#' and specifies that numeric variables from the clinical data table will be
#' used.
#' @param ... Variables not used by the \code{matrix} nor the
#' \code{MultiAssayExperiment} method which are passed into and used by the
#' \code{DataFrame} method.
#' @param verbose Default: 3. A number between 0 and 3 for the amount of
#' progress messages to give.  This function only prints progress messages if
#' the value is 3.
#' @return A vector of feature indicies, from the most promising features in the
#' first position to the least promising feature in the last position.
#' @importFrom survival coxph
#' @rdname coxphRanking
#' @usage NULL
#' @export
setGeneric("coxphRanking", function(measurementsTrain, ...) standardGeneric("coxphRanking"))

#' @rdname coxphRanking
#' @export
setMethod("coxphRanking", "matrix", function(measurementsTrain, survivalTrain, ...) # Matrix of numeric measurements.
{
  coxphRanking(S4Vectors::DataFrame(measurementsTrain, check.names = FALSE), survivalTrain, ...)
})

#' @rdname coxphRanking
#' @export
setMethod("coxphRanking", "DataFrame", function(measurementsTrain, survivalTrain, verbose = 3) # Clinical data or one of the other inputs, transformed.
{
  splitDataset <- .splitDataAndOutcome(measurementsTrain, survivalTrain)
  measurementsTrain <- splitDataset[["measurements"]]
  survivalTrain <- splitDataset[["outcome"]]

  if(any(sapply(measurementsTrain, class) %in% c("character", "factor"))){
  pValues <- apply(measurementsTrain, 2, function(featureColumn){
    fit <- survival::coxph(survivalTrain ~ featureColumn)
    s <- summary(fit)
    s$waldtest["pvalue"]
  })
  }else{
  tests <- colCoxTests(as.matrix(measurementsTrain), survivalTrain)
  pValues <- tests[colnames(measurementsTrain), "p.value"]
  }
  
  order(pValues) # From smallest to largest.
})

# One or more omics data sets, possibly with clinical data.
#' @rdname coxphRanking
#' @export
setMethod("coxphRanking", "MultiAssayExperiment", function(measurementsTrain, targets = names(measurementsTrain), survivalTrain, ...)
{
  tablesAndSurvival <- .MAEtoWideTable(measurementsTrain, targets, survivalTrain)
  measurementsTrain <- tablesAndSurvival[["dataTable"]]
  survivalTrain <- tablesAndSurvival[["outcome"]]
  
  if(ncol(measurementsTrain) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    coxphRanking(measurementsTrain, survivalTrain, ...)
})




# Copied verbatim from the archived survHD package Christoph Bernau and Levi Waldron
# Available on bitbucket https://bitbucket.org/lwaldron/survhd/src/master/R/filter.r
###############################################################################
# Filename: filter.r
# Title: Gene selection (filter) methods.
#
# Author: Christoph Bernau
# Email: <bernau@ibe.med.uni-muenchen.de>
# Date of creation: 28.8.2012
#
# Brief description:
#   Returns an object of class 'GeneSel'.
#
# Further comments and notes:
#   Are usually not called directly by the user, but via
#   'GeneSelection'.
#
###############################################################################

#' @useDynLib ClassifyR, .registration = TRUE
coxmatC<-function(X,time,status){
  ### method for handling ties (alternative 'breslow')
  method <- "efron"
  result<-	.C("coxmat", regmat = as.double(X), ncolmat = 
                as.integer(ncol(X)), nrowmat = as.integer(nrow(X)), 
              reg = as.double(X[, 1]), zscores = as.double(numeric(ncol(X))),
              coefs = as.double(numeric(ncol(X))), maxiter = as.integer(20), 
              nusedx = as.integer(nrow(X)), nvarx = as.integer(1), 
              time = as.double(time), status = as.integer(status), 
              offset = as.double(numeric(nrow(X))), 
              weights = as.double(numeric(nrow(X)) + 1), 
              strata = as.integer(numeric(nrow(X))), means = double(1), 
              beta = double(1), u = double(1), imat = double(1), 
              loglik = double(2), flag = integer(1), work = 
                double(2 * nrow(X) + 2 + 3), eps = as.double(1e-09), 
              tol_chol = as.double(.Machine$double.eps^0.75), 
              sctest = as.double(method == "efron"), sctest2 = as.double(1), 
              sctest3 = as.double(1), PACKAGE = "ClassifyR")	
  return(result)	
}


fastCox <- function(X, y, learnind, criterion, ...) {
  ### use learningset only and sort according to time
  X <- X[learnind, ]
  time <- y[learnind, 1]
  status <- y[learnind, 2]
  sorted <- order(time)
  time <- time[sorted]
  status <- status[sorted]
  X <- as.matrix(X[sorted, ])
  # compute columnwise coxmodels
  out <- coxmatC(X,time,status)
  # compute p-values
  if (criterion == "pvalue") 
    crit <- (1 - pnorm(abs(out$zscores))) * 2
  if (criterion == "coefficient") 
    crit <- abs(out$coefs)
  
  ### and return a VarSelOut-object
  new("VarSelOut", varsel = crit, criterion = criterion)
}

# equivalent to genefilter::rowttests for the cox model.  This is much faster
# than calling coxph for each row of a ##igh-dimensional matrix.
colCoxTests <- function(X, y, option = c("fast", "slow"), ...) {
  option <- match.arg(option)
  if (identical(option, "fast")) {
    X <- (as.matrix(X))  #make variables columns
    time <- y[, 1]
    status <- y[, 2]
    sorted <- order(time)
    time <- time[sorted]
    status <- status[sorted]
    X <- X[sorted, ]
    ## method for handling ties (alternative 'breslow')
    method <- "efron"
    ## compute columnwise coxmodels
    out <- coxmatC(X,time,status)
    ## compute p-values and return them
    output <- data.frame(coef = out$coefs, se.coef = out$coefs/out$zscores, 
                         p.value = (1 - 
                                      pnorm(abs(out$zscores))) * 2) 
    rownames(output) <- colnames(X)
  } else if (identical(option, "slow")) {
    output <- (apply(X, 2, function(xrow) {
      fit <- try(coxph(y ~ xrow))
      if (class(fit) == "try-error") {
        c(NA, NA)
      } else {
        summary(fit)$coefficients[1, c(1, 3, 5)]
      }
    }))
    colnames(output) <- c("coef", "se.coef", "p.value")
    rownames(output) <- rownames(X)
    output <- data.frame(output)
  } else stop("rowCoxTests: option should be fast or slow.")
  return(output)
  ### dataframe with two columns: coef = Cox regression
  ###coefficients, p.value =
  ### Wald Test p-values.  Rows correspond to the rows of X.
}