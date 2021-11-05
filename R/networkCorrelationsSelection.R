setGeneric("networkCorrelationsSelection", function(measurements, ...)
           standardGeneric("networkCorrelationsSelection"))

setMethod("networkCorrelationsSelection", "matrix", # Matrix of numeric measurements.
          function(measurements, classes, metaFeatures = NULL, ...)
{
  if(is.null(metaFeatures))
    stop("'metaFeatures' is NULL but must be provided.")

  networkCorrelationsSelection(DataFrame(t(measurements), check.names = FALSE), classes,
                               metaFeatures, ...)
})

# metaFeatures must also be a DataFrame.
setMethod("networkCorrelationsSelection", "DataFrame", # Possibly mixed data types.
          function(measurements, classes, metaFeatures = NULL, featureSets,
                   datasetName, trainParams, predictParams, resubstituteParams,
                   selectionName = "Differential Correlation of Sub-networks", verbose = 3)
{
  if(is.null(metaFeatures))
    stop("'metaFeatures' is NULL but must be provided.")
            
  if(verbose == 3)
    message("Ranking sub-networks by differences in correlation.")            

  networkIDs <- names(featureSets@sets)
  edgesPerNetwork <- sapply(featureSets@sets, nrow)
  networkIDsPerEdge <- factor(rep(networkIDs, edgesPerNetwork), levels = networkIDs)
  allInteractions <- do.call(rbind, featureSets@sets)
  
  if(!all(allInteractions[, 1] %in% colnames(measurements)) && !all(allInteractions[, 2] %in% colnames(measurements)))
    stop("Some interactors are not found in 'measurements'. Ensure that 'featureSets' does not have
       any features not in 'measurements'.")
  
  interactorTable <- measurements[, allInteractions[, 1]]
  otherInteractorTable <- measurements[, allInteractions[, 2]]
  
  oneClassTraining <- which(classes == levels(classes)[1])
  otherClassTraining <- which(classes == levels(classes)[2])
  oneClassFeatureCorrelations <- mapply(cor, interactorTable[oneClassTraining, ], otherInteractorTable[oneClassTraining, ])
  otherClassFeatureCorrelations <- mapply(cor, interactorTable[otherClassTraining, ], otherInteractorTable[otherClassTraining, ])
  oneClassSubnetworkCorrelations <- sapply(split(oneClassFeatureCorrelations, networkIDsPerEdge), mean)
  otherClassSubnetworkCorrelations <- sapply(split(otherClassFeatureCorrelations, networkIDsPerEdge), mean)
  overallSubnetworkCorrelations <- mapply(function(oneCorr, otherCorr) mean(c(oneCorr, otherCorr)), oneClassSubnetworkCorrelations, otherClassSubnetworkCorrelations)

  BSS <- mapply(function(totalEdges, oneClassCorr, otherClassCorr, overallCorr)
         {
           totalEdges * (oneClassCorr - overallCorr)^2 + totalEdges * (otherClassCorr - overallCorr)^2
         }, as.list(edgesPerNetwork), oneClassSubnetworkCorrelations, otherClassSubnetworkCorrelations, overallSubnetworkCorrelations)
  
  WSS <- mapply(function(oneClassFeatureCorr, otherClassFeatureCorr, oneClassNetworkCorr, otherClassNetworkCorr)
         {
           sum((oneClassFeatureCorr - oneClassNetworkCorr)^2) + sum((otherClassFeatureCorr - otherClassNetworkCorr)^2)
         }, split(oneClassFeatureCorrelations, networkIDsPerEdge), split(otherClassFeatureCorrelations, networkIDsPerEdge),
         as.list(oneClassSubnetworkCorrelations), as.list(otherClassSubnetworkCorrelations))

  networkRanking <- BSS/WSS
  orderedNetworks <- order(networkRanking, decreasing = TRUE)
  .pickFeatures(metaFeatures, classes, featureSets, trainParams, predictParams,
                resubstituteParams, orderedNetworks, verbose)            
})

setMethod("networkCorrelationsSelection", "MultiAssayExperiment", # Pick one numeric table from the data set.
          function(measurements, target = NULL, metaFeatures = NULL, ...)
{
  if(is.null(metaFeatures))
    stop("'metaFeatures' is NULL but must be provided.")            
            
  tablesAndClasses <- .MAEtoWideTable(measurements, target)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]            
  networkCorrelationsSelection(dataTable, classes, ...)
})