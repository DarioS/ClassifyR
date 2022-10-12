# Random Forest
RFparams <- function() {
    trainParams <- TrainParams(randomForestTrainInterface, tuneParams = list(mTryProportion = c(0.25, 0.33, 0.50, 0.66, 0.75, 1.00), ntree = seq(100, 500, 100)),
                               getFeatures = forestFeatures)
    predictParams <- PredictParams(randomForestPredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Random Survival Forest
RSFparams <- function() {
    trainParams <- TrainParams(rfsrcTrainInterface, tuneParams = list(mTryProportion = c(0.25, 0.33, 0.50, 0.66, 0.75, 1.00), ntree = seq(100, 500, 100)))
    predictParams <- PredictParams(rfsrcPredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

XGBparams <- function()
{
    trainParams <- TrainParams(extremeGradientBoostingTrainInterface, tuneParams = list(mTryProportion = c(0.25, 0.33, 0.50, 0.66, 0.75, 1.00), nrounds = c(5, 10, 15)),
                               getFeatures = XGBfeatures)
    predictParams <- PredictParams(extremeGradientBoostingPredictInterface)
    
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
    trainParams <- TrainParams(elasticNetGLMtrainInterface, getFeatures = elasticNetFeatures)
    predictParams <- PredictParams(elasticNetGLMpredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Support Vector Machine
SVMparams = function() {
    trainParams <- TrainParams(SVMtrainInterface,
                   tuneParams = list(kernel = c("linear", "polynomial", "radial", "sigmoid"),
                                     cost = 10^(-3:3)))
    predictParams <- PredictParams(SVMpredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Nearest Shrunken Centroid
NSCparams = function() {
    trainParams <- TrainParams(NSCtrainInterface, getFeatures = NSCfeatures)
    predictParams <- PredictParams(NSCpredictInterface)
    
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
    trainParams <- TrainParams(naiveBayesKernel, tuneParams = list(difference = c("unweighted", "weighted")))
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