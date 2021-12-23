# RandomForest
#' @importFrom ClassifyR TrainParams PredictParams randomForestTrainInterface randomForestPredictInterface forestFeatures
rfParams <- function() {
    trainParams = ClassifyR::TrainParams(
        ClassifyR::randomForestTrainInterface
    )
    
    predictParams = ClassifyR::PredictParams(predictor = ClassifyR::randomForestPredictInterface,
                                             returnType = "both")
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# ElatsicNet
#' @importFrom ClassifyR TrainParams PredictParams elasticNetGLMtrainInterface elasticNetGLMpredictInterface
elasticParams <- function() {
    trainParams = ClassifyR::TrainParams(ClassifyR::elasticNetGLMtrainInterface)
    
    predictParams = ClassifyR::PredictParams(predictor = ClassifyR::elasticNetGLMpredictInterface,
                                             returnType = "both")
    
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

# Logistic Params
#' @importFrom ClassifyR TrainParams PredictParams
logisticParams <- function() {
    trainParams = ClassifyR::TrainParams(logisticTrainInterface)
    predictParams = ClassifyR::PredictParams(predictor = logisticPredictInterface)
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

#' @importFrom ClassifyR TrainParams PredictParams SVMtrainInterface SVMpredictInterface
svmParams = function() {
    trainParams = ClassifyR::TrainParams(ClassifyR::SVMtrainInterface)
    predictParams = ClassifyR::PredictParams(predictor = ClassifyR::SVMpredictInterface,
                                             returnType = "both")
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

#' @importFrom ClassifyR TrainParams PredictParams DLDAtrainInterface DLDApredictInterface
DLDAParams = function() {
    trainParams = ClassifyR::TrainParams(ClassifyR::DLDAtrainInterface)
    predictParams = ClassifyR::PredictParams(predictor = ClassifyR::DLDApredictInterface,
                                             returnType = "both")
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}

#' @importFrom ClassifyR TrainParams PredictParams naiveBayesKernel
naiveBayesParams = function() {
    trainParams = ClassifyR::TrainParams(ClassifyR::naiveBayesKernel)
    predictParams = PredictParams(
        predictor = NULL,
        weighted = "weighted",
        weight = "height difference",
        returnType = "both"
    )
    return(list(trainParams = trainParams, predictParams = predictParams))
}

#' @importFrom ClassifyR TrainParams PredictParams elasticNetGLMpredictInterface
elasticNetPreval = function() {
    trainParams = ClassifyR::TrainParams(elasticNetGLMtrainInterfacePreval)
    
    predictParams = ClassifyR::PredictParams(predictor = ClassifyR::elasticNetGLMpredictInterface,
                                             returnType = "both")
    
    
    return(list(trainParams = trainParams, predictParams = predictParams))
}