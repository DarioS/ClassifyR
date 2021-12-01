setGeneric("NSCfeatures", function(trained, measurements, classes)
standardGeneric("NSCfeatures"))

setMethod("NSCfeatures", "pamrtrained",
          function(trained, measurements, classes)
{
  if(!requireNamespace("pamr", quietly = TRUE))
    stop("The package 'pamr' could not be found. Please install it.")
  
  minError <- min(trained[["errors"]])
  threshold <- trained[["threshold"]][max(which(trained[["errors"]] == minError))]
  params <- c(list(trained), list(list(x = t(as.matrix(measurements)), y = classes, geneid = 1:ncol(measurements))), threshold)
  chosen <- as.numeric(do.call(pamr::pamr.listgenes, params)[, 1])
  
  if(is.null(S4Vectors::mcols(measurements)))
    chosen <- colnames(measurements)[chosen]
  else
    chosen <- S4Vectors::mcols(measurements)[chosen, ]
  
  list(NULL, chosen)
})