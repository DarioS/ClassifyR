# RandomForest
RFparams <- function() {
    trainParams <- TrainParams(randomForestTrainInterface)
    predictParams <- PredictParams(randomForestPredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# k Nearest Neighbours
kNNparams <- function() {
    trainParams <- TrainParams(kNNinterface, tuneParams = list(k = 1:5))
    predictParams <- NULL
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Ordinary GLM
GLMparams <- function() {
    trainParams <- TrainParams(GLMtrainInterface)
    predictParams <- PredictParams(GLMpredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Elastic net GLM
elasticNetGLMparams <- function() {
    trainParams <- TrainParams(elasticNetGLMtrainInterface)
    predictParams <- PredictParams(elasticNetGLMpredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Support Vector Machine
SVMparams = function() {
    trainParams <- TrainParams(SVMtrainInterface)
    predictParams <- PredictParams(SVMpredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Diagonal Linear Discriminant Analysis
DLDAparams = function() {
    trainParams <- TrainParams(DLDAtrainInterface)
    predictParams <- PredictParams(DLDApredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# naive Bayes Kernel
naiveBayesParams <- function() {
    trainParams <- TrainParams(naiveBayesKernel)
    predictParams <- NULL
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Mixtures of Normals
mixModelsParams <- function() {
    trainParams <- TrainParams(mixModelsTrain, nbCluster = 1:2)
    predictParams <- PredictParams(mixModelsPredict)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Cox Proportional Hazards Model for Survival
coxphParams <- function() {
    trainParams <- TrainParams(coxphTrainInterface)
    predictParams <- PredictParams(predictor = coxphPredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Cox Proportional Hazards Model with Elastic Net for Survival
coxnetParams <- function() {
    trainParams <- TrainParams(coxnetTrainInterface)
    predictParams <- PredictParams(coxnetPredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}