# RandomForest
rfParams <- function() {
    trainParams = TrainParams(
        randomForestTrainInterface
    )
    
    predictParams = PredictParams(predictor = randomForestPredictInterface,
                                             returnType = "both")
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# ElatsicNet
elasticParams <- function() {
    trainParams = TrainParams(elasticNetGLMtrainInterface)
    
    predictParams = PredictParams(predictor = elasticNetGLMpredictInterface,
                                             returnType = "both")
    
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Logistic Params
logisticParams <- function() {
    trainParams = TrainParams(logisticTrainInterface)
    predictParams = PredictParams(predictor = logisticPredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

svmParams = function() {
    trainParams = TrainParams(SVMtrainInterface)
    predictParams = PredictParams(predictor = SVMpredictInterface,
                                             returnType = "both")
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

DLDAParams = function() {
    trainParams = TrainParams(DLDAtrainInterface)
    predictParams = PredictParams(predictor = DLDApredictInterface,
                                             returnType = "both")
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

naiveBayesParams = function() {
    trainParams = TrainParams(naiveBayesKernel)
    predictParams = PredictParams(
        predictor = NULL,
        weighted = "weighted",
        weight = "height difference",
        returnType = "both"
    )
    return(list(trainParams = trainParams, predictParams = predictParams))
}

mixModelsParams = function() {
    trainParams = TrainParams(mixModelsTrain)
    predictParams = PredictParams(mixModelsPredict)
    return(list(trainParams = trainParams, predictParams = predictParams))
}

elasticNetPreval = function() {
    trainParams = TrainParams(elasticNetGLMtrainInterfacePreval)
    
    predictParams = PredictParams(predictor = elasticNetGLMpredictInterface,
                                             returnType = "both")
    
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}


coxphParams = function() {
    trainParams = TrainParams(coxphTrainInterface)
    predictParams = PredictParams(predictor = coxphPredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}


coxnetParams = function() {
    trainParams = TrainParams(coxnetTrainInterface)
    predictParams = PredictParams(predictor = coxnetPredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

