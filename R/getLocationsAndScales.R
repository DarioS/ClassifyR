# Calculate Location and Scale

getLocationsAndScales <- function(measurements, location = c("mean", "median"), scale = c("SD", "MAD", "Qn"))
{
  location <- match.arg(location) 
  scale <- match.arg(scale)
   
  if(scale == "Qn" && !requireNamespace("robustbase", quietly = TRUE))
    stop("The package 'robustbase' could not be found. Please install it.")
 
  setNames(list(switch(location,
                       mean = apply(measurements, 2, mean),
                       median = apply(measurements, 2, median)),
                switch(scale,
                       SD = apply(measurements, 2, sd),
                       MAD = apply(measurements, 2, mad),
                       Qn = apply(measurements, 2, robustbase::Qn))),
                c(location, scale))
}