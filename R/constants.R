.ClassifyRenvir <- new.env(parent = emptyenv())

# Used internally during parameter selection based on best performance.
.ClassifyRenvir[["performanceInfoTable"]] <- matrix(c("Error", "lower",
                                                      "Accuracy", "higher",
                                                      "Balanced Error", "lower",
                                                      "Balanced Accuracy", "higher",
                                                      "Micro Precision", "higher",
                                                      "Micro Recall", "higher",
                                                      "Micro F1", "higher",
                                                      "Macro Precision", "higher",
                                                      "Macro Recall", "higher",
                                                      "Macro F1", "higher",
                                                      "Matthews Correlation Coefficient", "higher",
                                                      "C-index", "higher"),
                                                    ncol = 2, byrow = TRUE, dimnames = list(NULL, c("type", "better"))
) |> as.data.frame()

# Nice-looking names for feature selection and classification functions, to automatically use
# in a variety of performance plots.
.ClassifyRenvir[["functionsTable"]] <- matrix(
  c("bartlettRanking", "Bartlett Test",
    "classifyInterface", "Poisson LDA",
    "differentMeansRanking", "Difference in Means",
    "DLDAtrainInterface", "Diagonal LDA",
    "DLDApredictInterface", "Diagonal LDA",
    "DMDranking", "Differences of Medians and Deviations",
    "edgeRranking", "edgeR LRT",
    "elasticNetGLMtrainInterface", "Elastic Net GLM",
    "elasticNetGLMpredictInterface", "Elastic Net GLM",
    "fisherDiscriminant", "Fisher's LDA",
    "kNNinterface", "k Nearest Neighbours",
    "KolmogorovSmirnovRanking", "Kolmogorov-Smirnov Test",
    "kTSPclassifier", "k Top-Scoring Pairs",
    "KullbackLeiblerRanking", "Kullback-Leibler Divergence",
    "leveneRanking", "Levene Test",
    "likelihoodRatioRanking", "Likelihood Ratio Test (Normal)",
    "limmaRanking", "Moderated t-test",
    "mixModelsTrain", "Mixtures of Normals",
    "mixModelsPredict", "Mixtures of Normals",
    "naiveBayesKernel", "Naive Bayes Kernel",
    "NSCtrainInterface",  "Nearest Shrunken Centroids",
    "NSCpredictInterface", "Nearest Shrunken Centroids",
    "pairsDifferencesRanking", "Pairs Differences",
    "previousSelection", "Previous Selection", 
    "previousTrained", "Previous Trained",
    "randomForestTrainInterface", "Random Forest",
    "randomForestPredictInterface", "Random Forest",
    "subtractFromLocation", "Location Subtraction",
    "SVMtrainInterface", "Support Vector Machine",
    "SVMpredictInterface", "Support Vector Machine",
    "coxphTrainInterface", "Cox Proportional Hazards",
    "coxphPredictInterface", "Cox Proportional Hazards",
    "coxphRanking", "Cox Proportional Hazards",
    "coxnetTrainInterface", "Penalised Cox Proportional Hazards",
    "coxnetPredictInterface", "Penalised Cox Proportional Hazards",
    "NEMOEtrainInterface", "Nutrition-Ecotype Mixture of Experts"),
  ncol = 2, byrow = TRUE, dimnames = list(NULL, c("character", "name"))
) |> as.data.frame()