setGeneric("networkCorrelationsSelection", function(measurements, ...)
           {standardGeneric("networkCorrelationsSelection")})

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
          function(measurements, classes, metaFeatures = NULL, networkSets,
                   datasetName, trainParams, predictParams, resubstituteParams,
                   selectionName = "Differential Correlation of Sub-networks", verbose = 3)
{
  if(is.null(metaFeatures))
    stop("'metaFeatures' is NULL but must be provided.")
            
  if(verbose == 3)
    message("Ranking sub-networks by differences in correlation.")            

  networkIDs <- names(networkSets@sets)
  edgesPerNetwork <- sapply(networkSets@sets, nrow)
  networkIDsPerEdge <- rep(names(networkSets@sets), edgesPerNetwork)
  allInteractions <- do.call(rbind, networkSets@sets)
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
  
  WSS <- mapply(function(oneClassFeatureCorr, otherClassFeatureCorr, oneClassOtherFeatureCorr, otherClassOtherFeatureCorr)
         {
           sum((oneClassFeatureCorr - oneClassOtherFeatureCorr)^2) + sum((otherClassFeatureCorr - otherClassOtherFeatureCorr)^2)
         }, split(oneClassFeatureCorrelations, networkIDsPerEdge), split(otherClassFeatureCorrelations, networkIDsPerEdge),
         as.list(oneClassSubnetworkCorrelations), as.list(otherClassSubnetworkCorrelations))

  networkRanking <- BSS/WSS
  orderedNetworks <- order(networkRanking, decreasing = TRUE)
  .pickFeatures(metaFeatures, classes, networkSets, datasetName, trainParams, predictParams,
                resubstituteParams, orderedNetworks, selectionName, verbose)            
})

setMethod("networkCorrelationsSelection", "MultiAssayExperiment", # Pick one numeric table from the data set.
          function(measurements, target = NULL, ...)
{
  tablesAndClasses <- .MAEtoWideTable(measurements, target)
  dataTable <- tablesAndClasses[["dataTable"]]
  classes <- tablesAndClasses[["classes"]]            
  networkCorrelationsSelection(dataTable, classes, ...)
})