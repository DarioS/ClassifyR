# An Interface for class Package's knn Function

kNNinterface <- function(measurementsTrain, classesTrain, measurementsTest, ..., verbose = 3)
{
  # Ensure same ordering for both tables.
  measurementsTest <- measurementsTest[, colnames(measurementsTrain), drop = FALSE]
  
  if(!requireNamespace("class", quietly = TRUE))
    stop("The package 'class' could not be found. Please install it.")
  if(verbose == 3)
    message("Fitting k Nearest Neighbours classifier to data and predicting classes.")
  
  class::knn(as.matrix(measurementsTrain), as.matrix(measurementsTest), classesTrain, ...)
}
attr(kNNinterface, "name") <- "kNNinterface"
