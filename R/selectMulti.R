selectMulti <- function(measurementsTrain, classesTrain, params, verbose = 0)
          {
              assaysIndices <- lapply(unique(S4Vectors::mcols(measurementsTrain)[["assay"]]), function(assay) which(S4Vectors::mcols(measurementsTrain)[["assay"]] == assay))
              assayTrain <- lapply(assaysIndices, function(assayIndices) measurementsTrain[, assayIndices])
              
              featuresIndices <- mapply(.doSelection, 
                                         measurements = assayTrain,
                                         modellingParams = params,
                                         MoreArgs = list(outcomeTrain = classesTrain, 
                                                         crossValParams = CrossValParams(permutations = 1, folds = 5), ###### Where to get this from?
                                                         verbose = 0), SIMPLIFY = FALSE
                                        )

              unlist(mapply(function(allDataIndices, withinIndices) allDataIndices[withinIndices],
                     assaysIndices, lapply(featuresIndices, "[[", 2), SIMPLIFY = FALSE))
}
attr(selectMulti, "name") <- "Union Selection"
