.ClassifyRenvir <- new.env(parent = emptyenv())

<<<<<<< HEAD
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
              "SVMtrainInterface", "SVMpredictInterface", "coxphTrainInterface", "coxphPredictInterface", "coxphRanking", "coxnetTrainInterface", "coxnetPredictInterface", "NEMoEtrainInterface", "NEMoEpredictInterface"),
               name = c("Bartlett Test", "Poisson LDA", "Difference in Means", "Diagonal LDA", "Diagonal LDA",
              "Differences of Medians and Deviations", "Easy-Hard Classifier", "Easy-Hard Classifier", "edgeR LRT",
              "Elastic Net GLM", "Elastic Net GLM", "Fisher's LDA", "k Nearest Neighbours",
              "Kolmogorov-Smirnov Test", "k Top-Scoring Pairs", "Kullback-Leibler Divergence",
              "Levene Test", "Likelihood Ratio Test (Normal)", "Moderated t-test",
              "Mixture Modelling", "Mixture Modelling", "Naive Bayes Kernel", "Nearest Shrunken Centroids",
              "Nearest Shrunken Centroids", "Pairs Differences", "Previous Selection", "Previous Trained",
              "Random Forest", "Random Forest", "Location Subtraction", "Support Vector Machine",
              "Support Vector Machine", "Cox proportional hazards", "Cox proportional hazards", "Cox proportional hazards", "Penalised cox proportional hazards", "Penalised cox proportional hazards", "NEMoE", "NEMoE")
              )
=======
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
                            "coxnetPredictInterface", "Penalised Cox Proportional Hazards"),
                            ncol = 2, byrow = TRUE, dimnames = list(NULL, c("character", "name"))
              ) |> as.data.frame()
>>>>>>> a1042463bd17b00e24b9ac6790be4ddda6bdd41d
