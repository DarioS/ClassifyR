.ClassifyRenvir <- new.env(parent = emptyenv())

.ClassifyRenvir[["performanceInfoTable"]] <- data.frame(type = c("error", "accuracy", "balanced error", "balanced accuracy",
                                            "micro precision", "micro recall",
                                            "micro F1", "macro precision",
                                            "macro recall", "macro F1", "matthews"),
                                             better = c("lower", "higher", "lower", "higher",
                                              "higher", "higher", "higher", "higher",
                                              "higher", "higher", "higher")
                                            )
# Used internally during parameter selection based on best performance.

.ClassifyRenvir[["functionsTable"]] <- data.frame(
              character = c("bartlettSelection", "classifyInterface", "differentMeansSelection",
              "DLDAtrainInterface", "DLDApredictInterface", "DMDselection",
              "easyHardClassifierTrain", "easyHardClassifierPredict", "edgeRselection", "elasticNetGLMtrainInterface",
              "elasticNetGLMpredictInterface", "fisherDiscriminant", "kNNinterface",
              "KolmogorovSmirnovSelection", "kTSPclassifier", "KullbackLeiblerSelection",
              "leveneSelection", "likelihoodRatioSelection", "limmaSelection",
              "mixModelsTrain", "mixModelsPredict", "naiveBayesKernel", "NSCtrainInterface",
              "NSCpredictInterface", "pairsDifferencesSelection", "previousSelection", "previousTrained",
              "randomForestTrainInterface", "randomForestPredictInterface", "subtractFromLocation",
              "SVMtrainInterface", "SVMpredictInterface", "networkCorrelationsSelection"),
               name = c("Bartlett Test", "Poisson LDA", "Difference in Means", "Diagonal LDA", "Diagonal LDA",
              "Differences of Medians and Deviations", "Easy-Hard Classifier", "Easy-Hard Classifier", "edgeR LRT",
              "Elastic Net GLM", "Elastic Net GLM", "Fisher's LDA", "k Nearest Neighbours",
              "Kolmogorov-Smirnov Test", "k Top-Scoring Pairs", "Kullback-Leibler Divergence",
              "Levene Test", "Likelihood Ratio Test (Normal)", "Moderated t-test",
              "Mixture Modelling", "Mixture Modelling", "Naive Bayes Kernel", "Nearest Shrunken Centroids",
              "Nearest Shrunken Centroids", "Pairs Differences", "Previous Selection", "Previous Trained",
              "Random Forest", "Random Forest", "Location Subtraction", "Support Vector Machine",
              "Support Vector Machine", "Differential Edge Correlations")
              )