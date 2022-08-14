selectMulti <- function(measurementsTrain, classesTrain, params, verbose = 0)
          {
              assayTrain <- sapply(unique(mcols(measurementsTrain)[["assay"]]), function(assay) measurementsTrain[, mcols(measurementsTrain)[["assay"]] %in% assay], simplify = FALSE)
              featuresIndices <- mapply(.doSelection, 
                                         measurements = assayTrain,
                                         modellingParams = params,
                                         MoreArgs = list(outcomeTrain = classesTrain, 
                                                         crossValParams = CrossValParams(permutations = 1, folds = 5), ###### Where to get this from?
                                                         verbose = 0), SIMPLIFY = FALSE
                                        )
              
              unique(unlist(lapply(featuresIndices, "[[", 2)))
          }
