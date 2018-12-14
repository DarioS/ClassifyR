edgesToHubNetworks <- function(edges, minCardinality = 5)
{
  if(class(edges) == "matrix")
  {
    allFeatures <- unique(as.vector(edges))
  } else { # data.frame
    if(!all(sapply(edges, class) == "character"))
      stop("Both columns must be character type.")
    allFeatures <- unique(unlist(edges))
  }

  featuresByFirstEdge <- split(edges[, 2], factor(edges[, 1], levels = allFeatures))
  featuresBySecondEdge <- split(edges[, 1], factor(edges[, 2], levels = allFeatures))
  featuresByHub <- mapply(c, featuresByFirstEdge, featuresBySecondEdge, SIMPLIFY = FALSE)
  featuresByHub <- lapply(featuresByHub, unique)
  featuresByHub <- featuresByHub[sapply(featuresByHub, length) >= minCardinality]
  if(length(featuresByHub) == 0)
    stop("No features have at least ", minCardinality, " interactors.")
  FeatureSetCollection(featuresByHub)
}
