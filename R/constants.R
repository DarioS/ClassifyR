.ClassifyRenvir <- new.env(parent = emptyenv())

.ClassifyRenvir[["performanceInfoTable"]] <- data.frame(type = c("Error", "Accuracy", "Balanced Error", "Balanced Accuracy",
                                            "Micro Precision", "Micro Recall",
                                            "Micro F1", "Macro precision",
                                            "Macro Recall", "Macro F1", "Matthews Correlation Coefficient", "C index"),
                                             better = c("lower", "higher", "lower", "higher",
                                              "higher", "higher", "higher", "higher",
                                              "higher", "higher", "higher", "higher")
                                            )
# Used internally during parameter selection based on best performance.

.ClassifyRenvir[["functionsTable"]] <- data.frame(
              character = c("bartlettRanking", "classifyInterface", "differentMeansRanking",
              "DLDAtrainInterface", "DLDApredictInterface", "DMDranking",
              "easyHardClassifierTrain", "easyHardClassifierPredict", "edgeRranking", "elasticNetGLMtrainInterface",
              "elasticNetGLMpredictInterface", "fisherDiscriminant", "kNNinterface",
              "KolmogorovSmirnovRanking", "kTSPclassifier", "KullbackLeiblerRanking",
              "leveneRanking", "likelihoodRatioRanking", "limmaRanking",
              "mixModelsTrain", "mixModelsPredict", "naiveBayesKernel", "NSCtrainInterface",
              "NSCpredictInterface", "pairsDifferencesRanking", "previousSelection", "previousTrained",
              "randomForestTrainInterface", "randomForestPredictInterface", "subtractFromLocation",
              "SVMtrainInterface", "SVMpredictInterface", "coxphTrainInterface", "coxphPredictInterface", "coxphRanking"),
               name = c("Bartlett Test", "Poisson LDA", "Difference in Means", "Diagonal LDA", "Diagonal LDA",
              "Differences of Medians and Deviations", "Easy-Hard Classifier", "Easy-Hard Classifier", "edgeR LRT",
              "Elastic Net GLM", "Elastic Net GLM", "Fisher's LDA", "k Nearest Neighbours",
              "Kolmogorov-Smirnov Test", "k Top-Scoring Pairs", "Kullback-Leibler Divergence",
              "Levene Test", "Likelihood Ratio Test (Normal)", "Moderated t-test",
              "Mixture Modelling", "Mixture Modelling", "Naive Bayes Kernel", "Nearest Shrunken Centroids",
              "Nearest Shrunken Centroids", "Pairs Differences", "Previous Selection", "Previous Trained",
              "Random Forest", "Random Forest", "Location Subtraction", "Support Vector Machine",
              "Support Vector Machine", "Cox proportional hazards", "Cox proportional hazards", "Cox proportional hazards")
              )