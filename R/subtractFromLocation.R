# Subtract Numeric Feature Measurements from a Per-Class Location Within Cross-validation
# Useful transformation for differential variability classification.

subtractFromLocation <- function(measurementsTrain, measurementsTest, location = c("mean", "median"),
                                 absolute = TRUE, verbose = 3)
{
  isNumeric <- sapply(measurementsTrain, is.numeric)
  if(sum(isNumeric) == 0)
    stop("No features are numeric but at least one must be.")
  
  if(sum(isNumeric) != ncol(measurementsTrain) && verbose == 3)
    message("Some columns are not numeric. Only subtracting values from location for\n", 
            "the columns containing numeric data.")
  
  location <- match.arg(location)
  measurementsTrain <- measurementsTrain[, isNumeric]
  measurementsTest <- measurementsTest[, isNumeric]
  if(location == "mean")
    locations <- apply(measurementsTrain, 2, mean, na.rm = TRUE)
  else # median.
    locations <- apply(measurementsTrain, 2, median, na.rm = TRUE)
  
  transformedTrain <- S4Vectors::DataFrame(t(apply(measurementsTrain, 1, '-', locations)))
  transformedTest <- S4Vectors::DataFrame(t(apply(measurementsTest, 1, '-', locations)))
  if(absolute == TRUE)
  {
    transformedTrain <- S4Vectors::DataFrame(lapply(transformedTrain, abs))
    transformedTest <- S4Vectors::DataFrame(lapply(transformedTest, abs))
  }
  
  if(verbose == 3)
    message("Subtraction from ", location,
            {if(absolute == TRUE) " and absolute transformation"}, " completed.")
  
  list(transformedTrain, transformedTest)
}
attr(subtractFromLocation, "name") <- "subtractFromLocation"