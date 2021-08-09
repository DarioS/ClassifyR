setGeneric("easyHardFeatures", function(easyHardClassifier, ...)
standardGeneric("easyHardFeatures"))

setMethod("easyHardFeatures", "EasyHardClassifier",
          function(easyHardClassifier)
          {
            selectedFeatures <- data.frame(dataset = character(0), features = character(0))
            if(!is.null(easyHardClassifier@easyClassifier))
              selectedFeatures <- unique(rbind(selectedFeatures, data.frame(dataset = unname(easyHardClassifier@datasetIDs["easy"]),
                                                                            feature = unname(unlist(lapply(easyHardClassifier@easyClassifier, '[', "feature"))))))
            if(!is.null(easyHardClassifier@hardClassifier) && !is.null(easyHardClassifier@hardClassifier[["selected"]]))
              selectedFeatures <- rbind(selectedFeatures, data.frame(dataset = unname(easyHardClassifier@datasetIDs["hard"]),
                                                                     feature = easyHardClassifier@hardClassifier[["selected"]]))
              
            list(NULL, selectedFeatures)
          })