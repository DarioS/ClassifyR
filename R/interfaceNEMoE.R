# data("PD")
# 
# Microbiome = DataFrame(PD$data_list$Microbiome$Phylum)
# Response = PD$data_list$Response
# Nutrition = DataFrame(PD$data_list$Nutrition)
# 
# mcols(Microbiome)$dataset <- c("microbiome")
# mcols(Nutrition)$dataset <- c("nutrition")
# 
# combine_data <- cbind(Microbiome, Nutrition)
# 
# Microbiome_matrix <- data.matrix(Microbiome)
# Nutrition_matrix <- data.matrix(Nutrition)
# 
# training_model <- NEMoEtrainInterface(combine_data, Response)
# pred <- NEMoEpredictInterface(training_model, test = combine_data)
# 

setGeneric("NEMoEtrainInterface", function(measurements, ...) standardGeneric("NEMoEtrainInterface"))

setMethod("NEMoEtrainInterface", "matrix", function(measurements, classes, ...)
{
  NEMoEtrainInterface(DataFrame(t(measurements), check.names = FALSE), classes, ...)
})

setMethod("NEMoEtrainInterface", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurements, classes, ..., verbose = 3)
          {
            splitDataset <- .splitDataAndClasses(measurements, classes)
            measurements <- splitDataset[["measurements"]]
            isNumeric <- sapply(measurements, is.numeric)
            measurements <- measurements[, isNumeric, drop = FALSE]
            
            if(sum(isNumeric) == 0)
              stop("No features are numeric but at least one must be.")
            
            if(!requireNamespace("NEMoE", quietly = TRUE))
              stop("The package 'NEMoE' could not be found. Please install it.")
            
            assayTrain <- sapply(unique(mcols(measurements)[["dataset"]]), function(x) measurements[,mcols(measurements)[["dataset"]]%in%x], simplify = FALSE)
            
            trainedNEMoEModel <- NEMoE::fitNEMoE(NEMoE::NEMoE_buildFromList(as.matrix(assayTrain[[1]]), as.matrix(assayTrain[[2]]), classes))
            
            if(verbose == 3)
              message("NEMoE training completed.")
            
            trainedNEMoEModel  
          })

setMethod("NEMoEtrainInterface", "MultiAssayExperiment",
          function(measurements, targets = names(measurements), classes, ...)
          {
            tablesAndClasses <- .MAEtoWideTable(measurements, targets, classes)
            measurements <- tablesAndClasses[["dataTable"]]
            classes <- tablesAndClasses[["classes"]]

            if(ncol(measurements) == 0)
              stop("No variables in data tables specified by \'targets\' are numeric.")
            else
              NEMoEtrainIterface(measurements, classes, ...)
          })

################################################################################
#
# Predict interface
#
################################################################################


setGeneric("NEMoEpredictInterface", function(trained_model, test, ...) standardGeneric("NEMoEpredictInterface"))

setMethod("NEMoEpredictInterface", c("NEMoE", "matrix"), function(trained_model, test, ...)
{
  NEMoEpredictInterface(trained, DataFrame(t(test), check.names = FALSE), ...)
})

setMethod("NEMoEpredictInterface", c("NEMoE", "DFrame"), function(trained_model, test, returnType = c("both", "class", "score"), verbose = 3)
{
  
  if(!requireNamespace("NEMoE", quietly = TRUE))
    stop("The package 'NEMoE' could not be found. Please install it.")
  
  isNumeric <- sapply(test, is.numeric)
  test <- test[, isNumeric, drop = FALSE]
  returnType <- match.arg(returnType)
  
  assayTrain <- sapply(unique(mcols(test)[["dataset"]]), function(x) test[,mcols(test)[["dataset"]]%in%x], simplify = FALSE)
  
  predictions <- NEMoE::NEMoE_predict(NEMoE = trained_model, X_new = as.matrix(assayTrain[[1]]), Z_new = as.matrix(assayTrain[[2]]))
  
  if(verbose == 3)
    message("NEMoE predictions made.")
  
  factors <- factor(ifelse(predictions$output > 0.5, "Y", "N"))
  
  #predictions
  score_matrix <- cbind(predictions$output, 1 - predictions$output)
  rownames(score_matrix) <- row.names(test)
  colnames(score_matrix) <- c("Y", "N")
  
  switch(returnType, class = factors, # Factor vector.
         score = score_matrix, # Numeric matrix.
         both = data.frame(class = factors, score_matrix, check.names = FALSE))
})

setMethod("NEMoEpredictInterface", c("NEMoE", "MultiAssayExperiment"), function(trained_model, test, targets = names(test), ...)
{
  test <- .MAEtoWideTable(test, targets)[["dataTable"]] # Remove any classes, if present.

  if(ncol(test) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    NEMoEpredictInterface(trained_model, test, ...)
})

