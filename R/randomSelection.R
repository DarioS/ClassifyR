# Automated Selection of Previously Selected Features
randomSelection <- function(measurementsTrain, classesTrain, nFeatures, verbose = 3)
{
  if(verbose == 3)
    message("Choosing random features.")

  sample(ncol(measurementsTrain), nFeatures) # Return indices, not identifiers.
}
attr(randomSelection, "name") <- "randomSelection"