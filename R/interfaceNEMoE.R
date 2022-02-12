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


################################################################################
#
# Train Interface
#
################################################################################


setGeneric("NEMOEtrainInterface", function(measurementsTrain, ...) standardGeneric("NEMOEtrainInterface"))

setMethod("NEMOEtrainInterface", "matrix", function(measurementsTrain, classesTrain, ...)
{
  NEMOEtrainInterface(DataFrame(t(measurementsTrain), check.names = FALSE), classesTrain, ...)
})

setMethod("NEMOEtrainInterface", "DataFrame", # Clinical data or one of the other inputs, transformed.
          function(measurementsTrain, classesTrain, ..., verbose = 3)
          {
            splitDataset <- .splitDataAndOutcomes(measurementsTrain, classesTrain)
            measurementsTrain <- splitDataset[["measurements"]]
            isNumeric <- sapply(measurementsTrain, is.numeric)
            measurementsTrain <- measurementsTrain[, isNumeric, drop = FALSE]
            
            if(sum(isNumeric) == 0)
              stop("No features are numeric but at least one must be.")
            
            if(!requireNamespace("NEMoE", quietly = TRUE))
              stop("The package 'NEMoE' could not be found. Please install it.")
            
            assayTrain <- sapply(unique(mcols(measurementsTrain)[["dataset"]]), function(x) measurementsTrain[,mcols(measurementsTrain)[["dataset"]]%in%x], simplify = FALSE)
            
            trainedNEMoEModel <- NEMoE::fitNEMoE(NEMoE::NEMoE_buildFromList(as.matrix(assayTrain[[1]]), as.matrix(assayTrain[[2]]), classesTrain))
            
            if(verbose == 3)
              message("NEMoE training completed.")
            
            trainedNEMoEModel  
          })

setMethod("NEMOEtrainInterface", "MultiAssayExperiment",
          function(measurementsTrain, targets = names(measurementsTrain), classesTrain, ...)
          {
            tablesAndClasses <- .MAEtoWideTable(measurementsTrain, targets, classesTrain)
            measurementsTrain <- tablesAndClasses[["dataTable"]]
            classesTrain <- tablesAndClasses[["classes"]]
            
            if(ncol(measurementsTrain) == 0)
              stop("No variables in data tables specified by \'targets\' are numeric.")
            else
              NEMOEtrainInterface(measurementsTrain, classesTrain, ...)
          })


################################################################################
#
# Predict interface
#
################################################################################


setGeneric("NEMOEpredictInterface", function(model, measurementsTest, ...) standardGeneric("NEMOEpredictInterface"))

setMethod("NEMOEpredictInterface", c("NEMoE", "matrix"), function(model, measurementsTest, ...)
{
  NEMOEpredictInterface(trained, DataFrame(t(measurementsTest), check.names = FALSE), ...)
})

setMethod("NEMOEpredictInterface", c("NEMoE", "DFrame"), function(model, measurementsTest, returnType = c("both", "class", "score"), verbose = 3)
{
  
  if(!requireNamespace("NEMoE", quietly = TRUE))
    stop("The package 'NEMoE' could not be found. Please install it.")
  
  isNumeric <- sapply(measurementsTest, is.numeric)
  measurementsTest <- measurementsTest[, isNumeric, drop = FALSE]
  returnType <- match.arg(returnType)
  
  assayTrain <- sapply(unique(mcols(measurementsTest)[["dataset"]]), function(x) measurementsTest[,mcols(measurementsTest)[["dataset"]]%in%x], simplify = FALSE)
  
  predictions <- NEMoE::NEMoE_predict(NEMoE = model, X_new = as.matrix(assayTrain[[1]]), Z_new = as.matrix(assayTrain[[2]]))
  
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

setMethod("NEMOEpredictInterface", c("NEMoE", "MultiAssayExperiment"), function(model, measurementsTest, targets = names(test), ...)
{
  measurementsTest <- .MAEtoWideTable(measurementsTest, targets)[["dataTable"]] # Remove any classes, if present.
  
  if(ncol(measurementsTest) == 0)
    stop("No variables in data tables specified by \'targets\' are numeric.")
  else
    NEMOEpredictInterface(model, measurementsTest, ...)
})